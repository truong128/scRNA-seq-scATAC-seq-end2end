#!/usr/bin/env python3

from __future__ import annotations
import argparse
import gzip
import html
import json
import logging
import os
from pathlib import Path
import re
import sys
import warnings

import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.backends.backend_pdf import PdfPages
import numpy as np
import pandas as pd
from scipy import sparse
import scanpy as sc
import yaml

LOG = logging.getLogger('io_scRNA')
plt.rcParams.update({'font.family':'DejaVu Sans', 'pdf.fonttype':42,
                     'ps.fonttype':42, 'savefig.dpi':300,
                     'axes.spines.top':False, 'axes.spines.right':False,
                     'axes.titleweight':'semibold', 'axes.labelsize':10,
                     'figure.facecolor':'white'})


def mkdir(path):
    Path(path).mkdir(parents=True, exist_ok=True)
    return Path(path)


def savefig(fig, folder, name, small=False):
    folder = mkdir(folder)
    fig.savefig(folder / (name + '.pdf'), bbox_inches='tight', facecolor='white')
    fig.savefig(folder / (name + '.png'), dpi=300 if not small else 180,
                bbox_inches='tight', facecolor='white')
    plt.close(fig)


def stamp(ax, txt='Exploratory analysis; independent validation required'):
    ax.text(0.0, -0.14, txt, transform=ax.transAxes,
            fontsize=7, color='#666666', ha='left', va='top')


def status(registry, name, state, note):
    registry.append({'step': name, 'status': state, 'detail': str(note)})
    LOG.info('%s: %s - %s', name, state, note)


def qc_boxplot(df, feature, output, name, threshold=None):
    ids = list(dict.fromkeys(df['sample_id'].astype(str)))
    fig, ax = plt.subplots(figsize=(max(7, len(ids) * 1.3), 4.5))
    items = [df.loc[df['sample_id'] == sid, feature].dropna().values for sid in ids]
    ax.boxplot(items, labels=ids, patch_artist=True,
               boxprops={'facecolor':'#BBD9E8'}, showfliers=False)
    if threshold is not None: ax.axhline(threshold, color='#D55E00', ls='--', lw=1.2)
    ax.set(ylabel=feature.replace('_',' '), title=name)
    ax.tick_params(axis='x', rotation=40)
    stamp(ax)
    savefig(fig, output, name)


def check_h5(path):
    import h5py
    with h5py.File(path, 'r') as f:
        if 'matrix' not in f or 'data' not in f['matrix']:
            raise RuntimeError(f'Not a valid modern 10x sparse matrix H5: {path}')


def measure_ambient(raw_path, sample_id, outdir):
    raw = sc.read_10x_h5(str(raw_path), gex_only=True)
    raw.var_names_make_unique()
    total = np.asarray(raw.X.sum(axis=1)).ravel()
    bg = raw.X[total <= 20, :]
    if bg.shape[0] < 30 or bg.nnz == 0:
        return {'sample_id': sample_id, 'background_barcodes': int(bg.shape[0]),
                'ambient_profile': 'insufficient background'}, None
    abundance = np.asarray(bg.sum(axis=0)).ravel()
    inds = np.argsort(abundance)[::-1][:30]
    profile = pd.DataFrame({'sample_id': sample_id,
                            'gene': raw.var_names[inds], 'background_umi': abundance[inds]})
    profile.to_csv(outdir / f'{sample_id}_ambient_background_genes.csv', index=False)
    stats = {'sample_id': sample_id, 'background_barcodes': int(bg.shape[0]),
             'background_umis': float(abundance.sum()), 'ambient_profile': 'measured_not_corrected'}
    del raw
    return stats, profile


def load_sample(row, folder, cfg, registry):
    sid = str(row.sample_id)
    filtered = folder / f'{sid}_filtered.h5'
    raw = folder / f'{sid}_raw.h5'
    for p in (filtered, raw):
        if not p.exists(): raise FileNotFoundError(p)
        check_h5(p)
    ad = sc.read_10x_h5(str(filtered), gex_only=True)
    ad.var_names_make_unique()
    ad.obs_names_make_unique()
    ad.obs['sample_id'] = sid
    ad.obs['patient_id'] = str(row.patient_id)
    ad.obs['condition'] = str(row.condition)
    ad.obs['orig_barcode'] = ad.obs_names.copy().astype(str)
    ad.var['mt'] = ad.var_names.str.startswith('MT-')
    ad.var['ribo'] = ad.var_names.str.match(r'^RP[SL]')
    ad.var['hb'] = ad.var_names.str.match(r'^HB[ABDEG]')
    sc.pp.calculate_qc_metrics(ad, qc_vars=['mt','ribo','hb'], percent_top=None, log1p=False, inplace=True)
    qc = ad.obs[['sample_id','patient_id','condition','total_counts','n_genes_by_counts',
                 'pct_counts_mt','pct_counts_ribo','pct_counts_hb']].copy()
    min_gene, max_gene, max_mt = cfg['min_genes'], cfg['max_genes'], cfg['max_mt']
    good = (ad.obs['n_genes_by_counts'] >= min_gene) & \
           (ad.obs['n_genes_by_counts'] <= max_gene) & \
           (ad.obs['pct_counts_mt'] <= max_mt)
    ad.obs['pass_qc'] = good.to_numpy()
    if cfg['run_scrublet'] and int(good.sum()) >= 100:
        try:
            candidate = ad[good].copy()
            with warnings.catch_warnings():
                warnings.simplefilter('ignore')
                sc.pp.scrublet(candidate, expected_doublet_rate=0.06, random_state=1234)
            ad.obs['doublet_score'] = np.nan
            ad.obs['predicted_doublet'] = False
            ad.obs.loc[candidate.obs_names, 'doublet_score'] = candidate.obs['doublet_score'].values
            ad.obs.loc[candidate.obs_names, 'predicted_doublet'] = candidate.obs['predicted_doublet'].astype(bool).values
            good = good & (~ad.obs['predicted_doublet'])
            n_pred = int(ad.obs['predicted_doublet'].sum())
            sc_vals = candidate.obs['doublet_score'].to_numpy(dtype=float)
            expected_n = 0.06 * len(sc_vals)
            if n_pred < 0.1 * expected_n:

                status(registry, f'{sid}:doublet', 'LOW_CONFIDENCE',
                       f'Scrublet predicted only {n_pred} doublets vs ~{expected_n:.0f} expected at 6%; '
                       f'score distribution not bimodal (median {np.median(sc_vals):.3f}, '
                       f'p99 {np.quantile(sc_vals, .99):.3f}, max {sc_vals.max():.3f}). '
                       'Doublet removal effectively NOT applied; T/NK doublets cannot be excluded. '
                       'Use scDblFinder or a manual threshold.')
            else:
                status(registry, f'{sid}:doublet', 'COMPLETE',
                       f'Scrublet predicted {n_pred} doublets')
        except Exception as exc:

            ad.obs['doublet_score'] = np.nan
            ad.obs['predicted_doublet'] = False
            status(registry, f'{sid}:doublet', 'FAILED_OPTIONAL', repr(exc))
    else:
        ad.obs['doublet_score'] = np.nan
        ad.obs['predicted_doublet'] = False
        status(registry, f'{sid}:doublet', 'SKIPPED', 'Disabled or fewer than 100 QC-passing barcodes')
    ad.obs['pass_final'] = good.to_numpy()
    qc['pass_qc'] = ad.obs['pass_qc'].values
    qc['pass_final'] = ad.obs['pass_final'].values
    qc['doublet_score'] = ad.obs['doublet_score'].values
    qc['predicted_doublet'] = ad.obs['predicted_doublet'].values
    retained = int(good.sum())
    if retained < 30: raise RuntimeError(f'Only {retained} cells retained for {sid}, check QC thresholds.')
    sub = ad[good].copy()
    sub.layers['counts'] = sub.X.copy().astype(np.float32)
    return sub, qc, {'sample_id': sid, 'n_input': ad.n_obs, 'n_retained': retained,
                     'fraction_retained': round(retained/ad.n_obs, 4),
                     'median_genes': round(float(qc['n_genes_by_counts'].median()),1),
                     'median_umis': round(float(qc['total_counts'].median()),1)}


def signature_scores(ad, signatures, min_present=2, prefix='signature'):
    found = {}
    for name, genes in signatures.items():
        genes = [g for g in genes if g in ad.var_names]
        if len(genes) < min_present:
            continue
        X = ad[:, genes].X
        v = np.asarray(X.mean(axis=1)).ravel()
        mu, sd = np.mean(v), np.std(v)
        ad.obs[f'{prefix}_{name}'] = (v-mu)/max(sd,1e-8)
        found[name] = genes
    return found


def annotate(ad, marker_file, output, registry):
    df = pd.read_csv(marker_file)
    markers = {r.cell_type: str(r.markers).split(';') for r in df.itertuples()}

    available = signature_scores(ad, markers, 2, prefix='marker')
    if len(available) < 3:
        raise RuntimeError('Too few canonical markers found. Check gene nomenclature (gene symbols expected).')
    scored = [f'marker_{x}' for x in available]
    med = ad.obs.groupby('leiden', observed=True)[scored].median()
    rows = []
    assignments = {}
    for cluster, row in med.iterrows():
        vals = row.sort_values(ascending=False)
        first = str(vals.index[0]).replace('marker_', '')
        margin = float(vals.iloc[0] - vals.iloc[1]) if len(vals)>1 else 0
        score = float(vals.iloc[0])
        pred = first if score > 0.1 and margin > 0.08 else 'Unresolved'
        assignments[str(cluster)] = pred
        rows.append({'leiden':cluster, 'putative_cell_type':pred,
                     'highest_marker_score': score, 'score_margin': margin,
                     'top_label':first, 'status':'unvalidated_marker_annotation'})
    ad.obs['cell_type'] = ad.obs['leiden'].astype(str).map(assignments).astype('category')
    pd.DataFrame(rows).to_csv(output/'cluster_annotation_confidence.csv', index=False)
    med.to_csv(output/'celltype_marker_scores_by_cluster.csv')
    status(registry,'marker_based_annotation','EXPLORATORY',
           'Putative cluster labels from canonical markers; review manually and validate with reference atlas')
    return markers


def umap_plot(ad, col, output, name, maxcells=18000):
    if 'X_umap' not in ad.obsm: return
    fig, ax=plt.subplots(figsize=(8.2,6.5))
    values=ad.obs[col].astype(str)
    rng=np.random.default_rng(42)
    ix=np.arange(ad.n_obs)
    if len(ix)>maxcells: ix=rng.choice(ix, maxcells, replace=False)
    vv=values.iloc[ix].values
    groups=list(pd.unique(vv))
    palette=plt.get_cmap('tab20', max(1,len(groups)))
    for k,g in enumerate(groups):
        mask=vv==g
        ax.scatter(ad.obsm['X_umap'][ix[mask],0],ad.obsm['X_umap'][ix[mask],1],
                   s=1.8, alpha=0.65, rasterized=True, label=g, color=palette(k))
    ax.set(xlabel='UMAP 1', ylabel='UMAP 2', title=name)
    ax.legend(loc='upper left', bbox_to_anchor=(1.02,1), markerscale=4,
              frameon=False, fontsize=8, title=col)
    stamp(ax)
    savefig(fig, output, name)


def feature_umap(ad, gene, outdir):
    if gene not in ad.var_names: return
    arr=np.asarray(ad[:,gene].X.toarray()).ravel()
    fig,ax=plt.subplots(figsize=(7,5))
    points=ax.scatter(ad.obsm['X_umap'][:,0],ad.obsm['X_umap'][:,1],
                      s=2, c=arr, cmap='viridis', rasterized=True)
    fig.colorbar(points,ax=ax, label=f'log normalized {gene}', fraction=.04)
    ax.set(title=f'{gene} expression',xlabel='UMAP 1',ylabel='UMAP 2')
    stamp(ax)
    savefig(fig,outdir,f'gene_{gene}')


def pathway_plot(ad, names, outdir):
    valid=[f'signature_{x}' for x in names if f'signature_{x}' in ad.obs]
    if not valid:return
    data=ad.obs.groupby('cell_type',observed=True)[valid].mean()
    fig,ax=plt.subplots(figsize=(max(9,len(valid)*1.3),max(4,0.5*len(data)+2)))
    im=ax.imshow(data.values,aspect='auto',cmap='coolwarm',vmin=-1.5,vmax=1.5)
    ax.set_xticks(range(len(data.columns)),[x.replace('signature_','').replace('_',' ') for x in data.columns],rotation=40,ha='right')
    ax.set_yticks(range(len(data)),data.index)
    ax.set(title='Immune transcriptional module scores (not protein activity)')
    fig.colorbar(im,ax=ax,shrink=.65,label='Mean standardized score')
    savefig(fig,outdir,'immune_signature_heatmap')
    data.to_csv(outdir.parent/'tables'/'immune_signature_mean_by_celltype.csv')


def cluster_markers(ad, outdir, figdir, registry):
    if ad.obs['leiden'].nunique()<2:
        status(registry,'cluster_markers','SKIPPED','Only one cluster')
        return
    sc.tl.rank_genes_groups(ad,groupby='leiden',method='wilcoxon',n_genes=min(40,ad.n_vars),
                            use_raw=False,pts=True)
    df=sc.get.rank_genes_groups_df(ad,group=None)
    df.to_csv(outdir/'cluster_marker_rankings_exploratory.csv',index=False)
    if not len(df):return

    unique=[]
    for group, group_df in df.groupby('group'):
        for gene in group_df['names'].head(3):
            if gene not in unique and gene in ad.var_names:
                unique.append(gene)
    if not unique:return
    groups=sorted(ad.obs['leiden'].astype(str).unique(),key=lambda s:int(s) if s.isdigit() else s)
    expr=[]
    for group in groups:
        idx=np.flatnonzero(ad.obs['leiden'].astype(str).to_numpy()==group)
        expr.append(np.asarray(ad[idx,unique].X.mean(axis=0)).ravel())
    z=np.asarray(expr)
    z=(z-z.mean(axis=0))/(z.std(axis=0)+1e-6)
    fig,ax=plt.subplots(figsize=(max(9,.48*len(unique)),max(5,.5*len(groups))))
    im=ax.imshow(np.clip(z,-2,2),aspect='auto',cmap='RdBu_r',vmin=-2,vmax=2)
    ax.set_xticks(range(len(unique)),unique,rotation=65,ha='right',fontsize=8)
    ax.set_yticks(range(len(groups)),groups)
    ax.set(title='Exploratory cluster marker heatmap',ylabel='Leiden cluster')
    fig.colorbar(im,ax=ax,shrink=.6,label='Gene-wise standardized log expression')
    savefig(fig,figdir,'cluster_marker_heatmap')
    status(registry,'cluster_markers','COMPLETE','Wilcoxon markers for annotation; do not interpret as between-patient DE')


def compositions(ad, outdir, registry):
    tab=pd.crosstab(ad.obs['sample_id'],ad.obs['cell_type'])
    frac=tab.div(tab.sum(axis=1),axis=0)
    tab.to_csv(outdir/'cell_counts_per_sample_type.csv')
    frac.to_csv(outdir/'cell_proportions_per_sample_type.csv')
    fig, ax=plt.subplots(figsize=(9,4.4))
    x=np.arange(len(tab))
    bottom=np.zeros(len(tab))
    palette=plt.get_cmap('tab20',max(1,len(tab.columns)))
    for i,c in enumerate(frac.columns):
        vals=frac[c].values
        ax.bar(x,vals,bottom=bottom,label=str(c),width=.7,color=palette(i))
        bottom+=vals
    ax.set_xticks(x,tab.index)
    ax.set(ylabel='Proportion among QC-passing cells',ylim=(0,1.05),title='Tumor TME cellular composition')
    ax.legend(loc='upper left',bbox_to_anchor=(1.02,1),frameon=False,fontsize=8)
    stamp(ax,'Descriptive only; cell recovery differs from absolute tissue abundance')
    savefig(fig,outdir.parent/'figures','celltype_composition')
    status(registry,'cell_abundance','DESCRIPTIVE_ONLY',
           'Per-sample proportions only. No compositional statistical test was run.')


def pseudobulk(ad, outdir, registry):
    keys=ad.obs[['sample_id','patient_id','condition','cell_type']].astype(str)
    genes=ad.var_names.to_numpy()
    mat=ad.layers['counts'].tocsr()
    summed=[]
    rec=[]
    for group, indexes in keys.groupby(['sample_id','patient_id','condition','cell_type'],sort=True).indices.items():
        if len(indexes)<10:continue
        counts=np.asarray(mat[indexes,:].sum(axis=0)).ravel().astype(np.int64)
        summed.append(counts)
        rec.append(dict(zip(['sample_id','patient_id','condition','cell_type'],group),n_cells=len(indexes)))
    if not summed:
        status(registry,'pseudobulk','SKIPPED','No strata with >=10 cells')
        return
    meta=pd.DataFrame(rec)
    cols=[f'{r.sample_id}|{r.cell_type}' for r in meta.itertuples()]
    counts=pd.DataFrame(np.column_stack(summed), index=genes, columns=cols)
    counts.index.name='gene'
    counts.to_csv(outdir/'pseudobulk_gene_by_sample_celltype.csv.gz',compression='gzip')
    meta.insert(0,'pseudobulk_id',cols)
    meta.to_csv(outdir/'pseudobulk_metadata.csv',index=False)
    status(registry,'pseudobulk','COMPLETE','Raw UMI sums; suitable as input to edgeR/DESeq2 after design validation')
    by_group=meta[['patient_id','condition']].drop_duplicates()
    levels=by_group.groupby('condition').patient_id.nunique()
    if len(levels)<2 or (levels<3).any():
        status(registry,'between_condition_DE','NOT_TESTABLE',
               f'Cannot make replicated between-group inference: donors by condition {levels.to_dict()}')
        (outdir/'no_DE_explanation.txt').write_text(
            'No between-condition DEG test. Need >=2 conditions and preferably >=3 independent biological donors per condition.\n'
            f'Observed independent donors by condition: {levels.to_dict()}.\n'
            'Do not compare individual cells as if they were independent patients.\n')
    else:
        status(registry,'between_condition_DE','REQUIRES_EDGER_MODULE',
               f'Donors per condition {levels.to_dict()}. The edgeR script can run a donor-level '
               'two-group test. Cell-type labels are marker-based and unvalidated, so those FDRs '
               'are not confirmed cell-type biology.')


def mechanism_population(label):
    text=str(label)
    return text != 'Unresolved' and 'TCRpos_ambiguous' not in text


def ligand_receptor(ad, pair_csv, outdir, registry):
    pairs=pd.read_csv(pair_csv)
    gene_index={x:i for i,x in enumerate(ad.var_names)}
    hits=[]

    for sample in sorted(ad.obs['sample_id'].astype(str).unique()):
        sub_obs=ad.obs[ad.obs['sample_id'].astype(str)==sample]
        counts=sub_obs['cell_type'].astype(str).value_counts()
        groups=[x for x in counts.index if counts.loc[x]>=30 and mechanism_population(x)]
        if len(groups)<2:continue
        means={}
        fractions={}
        for group in groups:
            ix=np.flatnonzero(((ad.obs['sample_id'].astype(str)==sample)&
                               (ad.obs['cell_type'].astype(str)==group)).to_numpy())
            X=ad.X[ix,:].tocsr()
            means[group]=np.asarray(X.mean(axis=0)).ravel()
            fractions[group]=np.asarray((X>0).mean(axis=0)).ravel()
        for p in pairs.itertuples():
            if p.ligand not in gene_index or p.receptor not in gene_index:continue
            i,j=gene_index[p.ligand],gene_index[p.receptor]
            for sender in groups:
                if fractions[sender][i]<.10:continue
                for receiver in groups:
                    if fractions[receiver][j]<.10:continue
                    score=np.sqrt(means[sender][i]*means[receiver][j])
                    hits.append({'sample_id':sample,'sender':sender,'receiver':receiver,
                                 'ligand':p.ligand,'receptor':p.receptor,'axis':p.axis,
                                 'ligand_fraction':fractions[sender][i],
                                 'receptor_fraction':fractions[receiver][j],
                                 'coexpression_score':float(score),
                                 'interpretation':'candidate_only_not_cellchat_or_causal'})
    hits=pd.DataFrame(hits)
    if not len(hits):
        status(registry,'cell_communication','NOT_TESTABLE',
               'No sample had >=2 confidently annotated populations (>=30 cells) with expressed pairs')
        return
    if len(hits):
        key=['sender','receiver','ligand','receptor','axis']
        hits['n_samples_supporting']=hits.groupby(key)['sample_id'].transform('nunique')
        hits=hits.sort_values(['n_samples_supporting','coexpression_score'],ascending=False)
        hits.to_csv(outdir/'ligand_receptor_candidate_coexpression.csv',index=False)
        top=hits.head(20).iloc[::-1]
        fig,ax=plt.subplots(figsize=(10,max(5,len(top)*.34)))
        labels=[f'{r.sample_id}: {r.sender} -> {r.receiver}: {r.axis}' for r in top.itertuples()]
        ax.barh(range(len(top)),top.coexpression_score.values,color='#5D8FA6')
        ax.set_yticks(range(len(top)),labels,fontsize=7)
        ax.set(xlabel='Heuristic co-expression score (NOT interaction probability)',
               title='Candidate ligand-receptor pairs')
        stamp(ax,'Hypotheses from expression only; no physical contact or signaling demonstrated')
        savefig(fig,outdir.parent/'figures','ligand_receptor_hypotheses')
    status(registry,'cell_communication','EXPLORATORY',
           f'{len(hits)} candidate sample-level pairs (computed per sample, not pooled across donors; '
           'see n_samples_supporting). Unresolved and TCRpos_ambiguous labels are excluded before pairing. '
           'Expression heuristic, not CellChat/LIANA/NicheNet inference')


def tcr_guard(ad, outdir, registry, min_tcr_fraction=0.20):
    t_labels = {'T_cells', 'CD8_T', 'CD4_T', 'Treg'}
    if 'tcr_clonotype' not in ad.obs or not (ad.obs['tcr_clonotype'].astype(str) != '').any():
        status(registry, 'tcr_guided_annotation', 'SKIPPED', 'No TCR data available to cross-check labels')
        return
    has_tcr = (ad.obs['tcr_clonotype'].astype(str) != '').to_numpy()
    cell_type = ad.obs['cell_type'].astype(str).to_numpy()
    leiden = ad.obs['leiden'].astype(str).to_numpy()
    detect = {}
    for gene in ['CD3D', 'CD3E', 'NCAM1']:
        if gene in ad.var_names:
            col = ad[:, gene].X
            col = col.toarray().ravel() if hasattr(col, 'toarray') else np.asarray(col).ravel()
            detect[gene] = col > 0
    rows, relabel = [], {}
    for cluster in sorted(set(leiden), key=lambda x: int(x) if x.isdigit() else x):
        m = leiden == cluster
        label = cell_type[m][0]
        frac = float(has_tcr[m].mean())
        row = {'leiden': cluster, 'marker_label': label, 'n_cells': int(m.sum()),
               'n_cells_with_tcr': int(has_tcr[m].sum()), 'tcr_fraction': round(frac, 4),
               'samples': ';'.join(f'{k}={v}' for k, v in
                                   ad.obs.loc[m, 'sample_id'].astype(str).value_counts().items())}
        mt = m & has_tcr
        for gene, flag in detect.items():
            row[f'{gene}_detected_in_tcr_cells'] = round(float(flag[mt].mean()), 4) if mt.any() else np.nan
        flagged = label not in t_labels and label != 'Unresolved' and frac >= min_tcr_fraction
        row['action'] = 'relabelled_ambiguous' if flagged else 'kept'
        if flagged:
            relabel[cluster] = f'{label}_TCRpos_ambiguous'
        rows.append(row)
    pd.DataFrame(rows).to_csv(outdir / 'tcr_guided_annotation_check.csv', index=False)
    if relabel:
        new = np.array([relabel.get(c, t) for c, t in zip(leiden, cell_type)])
        ad.obs['cell_type_marker_only'] = ad.obs['cell_type'].astype(str).to_numpy()
        ad.obs['cell_type'] = pd.Categorical(new)
        n = int(sum(np.isin(leiden, list(relabel))))
        status(registry, 'tcr_guided_annotation', 'RELABELLED',
               f'{len(relabel)} cluster(s) ({n} cells) labelled non-T by markers but >= {min_tcr_fraction:.0%} '
               f'TCR+: {relabel}. Need reference/CITE-seq/doublet review; see tcr_guided_annotation_check.csv')
    else:
        status(registry, 'tcr_guided_annotation', 'COMPLETE',
               f'No non-T cluster has >= {min_tcr_fraction:.0%} TCR+ cells')


def clonotypes(ad, inp, samples, outdir, registry):
    summary=[]
    cc=[]
    ad.obs['tcr_clonotype']=''
    for row in samples.itertuples(index=False):
        sid=str(row.sample_id)
        tcr=inp/f'{sid}_tcr.csv.gz'
        if not tcr.exists(): continue
        with gzip.open(tcr,'rt') as f:
            df=pd.read_csv(f)
        if len(df)==0 or 'raw_clonotype_id' not in df or 'barcode' not in df:
            status(registry,f'{sid}:tcr','SKIPPED','No actual TCR contigs available')
            continue
        if 'productive' in df:
            df=df[df['productive'].astype(str).str.lower()=='true'].copy()
        if 'chain' in df:
            df=df[df['chain'].astype(str).str.startswith('TR')].copy()
        df=df.dropna(subset=['raw_clonotype_id'])

        df=df[['barcode','raw_clonotype_id']].drop_duplicates('barcode')
        barcodes=ad.obs.loc[ad.obs['sample_id']==sid,'orig_barcode'].astype(str)
        idmap=df.set_index('barcode')['raw_clonotype_id'].astype(str)
        matched=barcodes.map(idmap)
        ad.obs.loc[barcodes.index,'tcr_clonotype']=matched.fillna('').to_numpy()
        observed=matched.dropna().value_counts()
        n_tcr=int(matched.notna().sum())
        summary.append({'sample_id':sid,'qc_cells':int(len(barcodes)),
                        'cells_with_tcr':n_tcr,'distinct_clonotypes':len(observed),
                        'largest_clone_size':int(observed.iloc[0]) if len(observed) else 0})
        for clone,n in observed.items(): cc.append({'sample_id':sid,'clonotype':clone,'cells':int(n)})
        status(registry,f'{sid}:tcr','COMPLETE',f'{n_tcr} matched QC cells with TCR contigs')
    if not summary:
        status(registry,'TCR_clonotypes','SKIPPED','No matched 5-prime TCR library data')
        return
    pd.DataFrame(summary).to_csv(outdir/'tcr_summary.csv',index=False)
    pd.DataFrame(cc).to_csv(outdir/'tcr_clonotype_counts.csv',index=False)
    fig,ax=plt.subplots(figsize=(6,4))
    x=[s['sample_id'] for s in summary]
    val=[s['cells_with_tcr'] for s in summary]
    ax.bar(x,val,color='#5A8D73')
    ax.set(ylabel='QC-passing cells with matched TCR',title='Matched TCR recovery')
    stamp(ax,'TCR expansion does not establish tumor antigen specificity')
    savefig(fig,outdir.parent/'figures','tcr_clonotype_recovery')
    status(registry,'TCR_clonotypes','COMPLETE','Joined TCR annotations by barcode and sample; clonal frequencies')


def trajectory(ad, output, registry):
    if 'TCF7' not in ad.var_names:
        status(registry,'T_cell_pseudotime','SKIPPED','TCF7 absent')
        return
    ix=(ad.obs['cell_type'].astype(str).str.contains('T_cells|CD8_T|CD4_T|Treg',regex=True)).to_numpy()
    if ix.sum()<100:
        status(registry,'T_cell_pseudotime','SKIPPED','Fewer than 100 putative T cells')
        return
    sub=ad[ix].copy()
    donor_frac=sub.obs['patient_id'].astype(str).value_counts(normalize=True)
    single_donor=float(donor_frac.iloc[0])>0.90
    donor_note=', '.join(f'{k}={v:.1%}' for k,v in donor_frac.items())
    try:
        nc=min(20,int(min(sub.n_vars,sub.n_obs)-1))
        sc.pp.pca(sub,n_comps=nc)
        sc.pp.neighbors(sub,n_neighbors=15,n_pcs=nc)
        sc.tl.diffmap(sub)
        root_vals=np.asarray(sub[:, 'TCF7'].X.toarray()).ravel()
        sub.uns['iroot']=int(np.argmax(root_vals))
        sc.tl.dpt(sub,n_dcs=min(10,sub.obsm['X_diffmap'].shape[1]-1))
        sub.obs[['sample_id','cell_type','dpt_pseudotime']].to_csv(output/'putative_Tcell_diffusion_pseudotime.csv')
        fig,ax=plt.subplots(figsize=(7,5))
        points=ax.scatter(sub.obsm['X_umap'][:,0],sub.obsm['X_umap'][:,1],
                          c=sub.obs['dpt_pseudotime'],s=3,cmap='plasma',rasterized=True)
        fig.colorbar(points,ax=ax,label='Diffusion pseudotime')
        ax.set(title='Exploratory T-cell state ordering',xlabel='UMAP 1',ylabel='UMAP 2')
        stamp(ax,'Computational ordering is not lineage or actual time evidence')
        savefig(fig,output.parent/'figures','Tcell_diffusion_pseudotime')
        if single_donor:
            status(registry,'T_cell_pseudotime','EXPLORATORY_SINGLE_DONOR',
                   f'Rooted at highest TCF7 cell. T cells are dominated by one donor ({donor_note}); '
                   'trajectory reflects that donor only and is not validated across patients')
        else:
            status(registry,'T_cell_pseudotime','EXPLORATORY',
                   f'Rooted at highest detected TCF7 cell. Donor composition: {donor_note}')
    except Exception as exc:
        status(registry,'T_cell_pseudotime','FAILED_OPTIONAL',repr(exc))


def immune_subclustering(ad, outdir, figdir, registry):
    groups={
        'T_lineage':['T_cells','CD8_T','CD4_T','Treg'],
        'Myeloid':['Monocytes','Macrophages','DC','pDC'],
        'B_lineage':['B_cells','Plasma'],
    }
    for label, celltypes in groups.items():
        col=f'subcluster_{label}'
        ad.obs[col]='not_in_lineage'
        index=np.flatnonzero(ad.obs['cell_type'].astype(str).isin(celltypes).to_numpy())
        if len(index)<100:
            status(registry,f'subcluster:{label}','SKIPPED',f'Only {len(index)} putative cells')
            continue
        try:
            subset=ad[index].copy()
            sc.pp.highly_variable_genes(subset,n_top_genes=min(1800,subset.n_vars-1),flavor='seurat')
            k=min(25,int(subset.var['highly_variable'].sum())-1,len(index)-2)
            sc.pp.pca(subset,n_comps=k,mask_var='highly_variable',random_state=1234)
            sc.pp.neighbors(subset,n_neighbors=12,n_pcs=k,random_state=1234)
            sc.tl.leiden(subset,resolution=.55,flavor='igraph',directed=False,random_state=1234)
            sc.tl.umap(subset,random_state=1234)
            ad.obs.loc[subset.obs_names,col]=[f'{label}_{i}' for i in subset.obs['leiden'].astype(str)]
            subset.obs[['sample_id','leiden','cell_type']].to_csv(outdir/f'{label}_subcluster_assignments.csv')
            umap_plot(subset,'leiden',figdir,f'{label}_subclustering')
            status(registry,f'subcluster:{label}','EXPLORATORY',
                   f'{len(index)} cells clustered; subclusters not yet biologically annotated')
        except Exception as exc:
            status(registry,f'subcluster:{label}','FAILED_OPTIONAL',repr(exc))


def signature_umap(ad, module, outdir):
    field=f'signature_{module}'
    if field not in ad.obs:return
    vals=ad.obs[field].to_numpy()
    fig,ax=plt.subplots(figsize=(7.0,5.3))
    p=ax.scatter(ad.obsm['X_umap'][:,0],ad.obsm['X_umap'][:,1],
                 c=vals,cmap='coolwarm',vmin=-2,vmax=2,s=2,rasterized=True)
    fig.colorbar(p,ax=ax,fraction=.04,label='Standardized expression module score')
    ax.set(xlabel='UMAP 1',ylabel='UMAP 2',title=module.replace('_',' '))
    stamp(ax,'RNA marker score, not a direct measure of signaling or effector function')
    savefig(fig,outdir,f'module_{module}')


def targets(ad,outdir,registry):
    genes=['CD274','PDCD1LG2','PVR','TIGIT','CSF1R','CSF1','CXCR4','CXCL12',
           'SPP1','TREM2','HAVCR2','LAG3','CTLA4','ENTPD1','TGFB1','B2M','HLA-A',
           'HLA-B','TAP1','TAP2','JAK1','JAK2','STAT3','IDO1','PDCD1','CD47']
    rows=[]
    for g in genes:
        if g not in ad.var_names:continue
        x=np.asarray(ad[:,g].X.toarray()).ravel()
        for typ, ids in ad.obs.groupby('cell_type',observed=True).indices.items():
            if not mechanism_population(typ):
                continue
            a=x[ids]
            rows.append({'gene':g,'putative_cell_type':str(typ),
                         'mean_log1p_expression':float(a.mean()),
                         'fraction_detected':float((a>0).mean()),
                         'evidence':'single_cell_expression_only_not_druggability'})
    pd.DataFrame(rows).to_csv(outdir/'candidate_target_expression_by_celltype.csv',index=False)
    status(registry,'target_prioritization','EXPLORATORY',
           'Expression-only candidate table. Unresolved and TCRpos_ambiguous labels are excluded. '
           'No causal, druggability, safety or therapeutic validation')


def write_report_html(outdir, registry, samples, n_obs, n_vars):
    donor_counts=samples.groupby(samples['condition'].astype(str))['sample_id'].nunique().sort_index()
    donor_txt=', '.join(f'{k} n={int(v)}' for k,v in donor_counts.items())
    table=''.join(f'<tr><td>{html.escape(str(r["step"]))}</td><td>{html.escape(str(r["status"]))}</td><td>{html.escape(str(r["detail"]))}</td></tr>' for r in registry)
    figlinks=''.join(f'<li><a href="figures/{html.escape(p.name)}">{html.escape(p.name)}</a></li>' for p in sorted((outdir/'figures').glob('*.pdf')))
    document=f'''<!DOCTYPE html><html lang="en"><head><meta charset="utf-8"><title>GSE274934 RNA, TCR, and ATAC report</title><style>
    body{{font-family:Arial,sans-serif;max-width:1160px;margin:30px auto;padding:0 24px;color:#14253d}}
    h1{{color:#174A69}}h2{{margin-top:2em;color:#174A69}}table{{border-collapse:collapse;font-size:13px;width:100%}}
    th,td{{padding:8px;text-align:left;border-bottom:1px solid #dadfe7;vertical-align:top}}
    th{{background:#eaf2f6}}.alert{{padding:16px;background:#fff1d1;border-left:5px solid #D18916}}
    </style></head><body><h1>GSE274934 RNA, TCR, and ATAC</h1>
    <div class="alert"><b>Scientific limitations:</b> Exploratory analysis. RNA donors by tumor type: {html.escape(donor_txt)}. Cell-type labels are marker-based and unvalidated. Where doublet status is LOW_CONFIDENCE, Scrublet did not remove the expected doublet rate. RNA and ATAC are separate libraries, so this report does not match an RNA cell to an ATAC cell and does not run WNN, MultiVI, or same-cell peak-to-gene linking. ATAC differential accessibility is not tested unless every tumor type has at least 3 ATAC donors. RNA edgeR, when the design has at least 3 donors in each tumor type, tests donors. Those FDRs are not confirmed cell-type biology. No batch correction is applied. There is no response annotation and no causal claim.</div>
    <h2>Dataset</h2><p>GEO GSE274934. RNA samples: {html.escape(', '.join(samples.sample_id.astype(str)))}. QC-retained RNA cells: {n_obs:,}. Genes after RNA filtering: {n_vars:,}. ATAC is analysed only for samples that have both a peak matrix and a Cell Ranger singlecell file. AL06 has RNA and TCR and no ATAC library.</p>
    <h2>Downloaded sample source</h2><p><a href="https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE274934">NCBI GEO GSE274934</a>. Tumor-type labels are the GEO sample characteristic <code>tumor type</code>, copied into the sample sheet.</p>
    <h2>Analysis and gating audit</h2><table><thead><tr><th>Module</th><th>Status</th><th>Reason</th></tr></thead><tbody>{table}</tbody></table>
    <h2>Vector figures (PDF)</h2><ul>{figlinks}</ul><h2>Machine-readable outputs</h2><p>tables/*.csv, tables/*.csv.gz, annotated_scrna.h5ad, cell_metadata.csv.gz, qc/*.csv, provenance.json.</p>
    <p><b>Before any manuscript claim:</b> confirm labels with a reference atlas, apply a validated doublet and ambient model, treat marker scores and edgeR FDRs as hypothesis-generating, and check an independent cohort plus an orthogonal assay.</p></body></html>'''
    (outdir/'index.html').write_text(document,encoding='utf-8')


def report(outdir, registry, samples, overview, ad):
    pd.DataFrame(registry).to_csv(outdir/'step_completeness_audit.csv',index=False)
    with (outdir/'step_completeness_audit.json').open('w') as f:
        json.dump(registry,f,indent=2)
    overview.to_csv(outdir/'sample_qc_summary.csv',index=False)
    ad.obs.to_csv(outdir/'cell_metadata.csv.gz',compression='gzip')
    ad.write_h5ad(outdir/'annotated_scrna.h5ad',compression='gzip')
    write_report_html(outdir, registry, samples, ad.n_obs, ad.n_vars)


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--samples',required=True)
    p.add_argument('--matrix-dir',required=True)
    p.add_argument('--config',required=True)
    p.add_argument('--markers',required=True)
    p.add_argument('--ligand-receptor',required=True)
    p.add_argument('--outdir',default='analysis')
    p.add_argument('--min-genes',type=int,default=200)
    p.add_argument('--max-genes',type=int,default=9000)
    p.add_argument('--max-mt',type=float,default=25)
    p.add_argument('--n-pcs',type=int,default=30)
    p.add_argument('--resolution',type=float,default=.6)
    p.add_argument('--run-scrublet',choices=['true','false'],default='true')
    p.add_argument('--gene-annotation',default='',
                   help='GENCODE GTF(.gz). Required for ATAC gene activity and RNA linking.')
    args=p.parse_args()
    log=logging.StreamHandler(sys.stdout)
    log.setFormatter(logging.Formatter('%(asctime)s %(levelname)s %(message)s'))
    LOG.setLevel(logging.INFO)
    LOG.addHandler(log)
    out=mkdir(args.outdir)
    figures=mkdir(out/'figures'); tables=mkdir(out/'tables'); qc_out=mkdir(out/'qc')
    inp=Path(args.matrix_dir)
    cfg={'min_genes':args.min_genes,'max_genes':args.max_genes,
         'max_mt':args.max_mt,'run_scrublet':args.run_scrublet=='true'}
    registry=[]
    with open(args.config) as f: analysis_cfg=yaml.safe_load(f)
    samples=pd.read_csv(args.samples,dtype=str).fillna('')
    if samples['sample_id'].duplicated().any(): raise ValueError('Duplicate sample_id in sample sheet')
    if samples['patient_id'].eq('').any(): raise ValueError('Missing patient_id')
    objects=[]; qcs=[]; overview=[]; ambient=[]
    for row in samples.itertuples(index=False):
        sid=str(row.sample_id)
        LOG.info('Loading %s',sid)
        a,q,ov=load_sample(row,inp,cfg,registry)
        objects.append(a)
        qcs.append(q)
        overview.append(ov)
        try:
            amb,_=measure_ambient(inp/f'{sid}_raw.h5',sid,qc_out)
            ambient.append(amb)
            status(registry,f'{sid}:ambient','MEASURED_ONLY','Unfiltered droplet profile, no subtraction performed')
        except Exception as exc:
            status(registry,f'{sid}:ambient','FAILED_OPTIONAL',repr(exc))
    all_qc=pd.concat(qcs,ignore_index=True)
    all_qc.to_csv(qc_out/'qc_all_cells.csv.gz',index=False,compression='gzip')
    pd.DataFrame(ambient).to_csv(qc_out/'ambient_assessment.csv',index=False)
    for key,name,limit in [('n_genes_by_counts','QC_gene_counts',args.min_genes),
                            ('total_counts','QC_UMIs',None),
                            ('pct_counts_mt','QC_mitochondrial_percent',args.max_mt)]:
        qc_boxplot(all_qc,key,figures,name,limit)
    status(registry,'QC','COMPLETE','Per-sample QC and transparent retention metrics')

    ad=sc.concat(objects,label=None,join='inner',index_unique='-')
    del objects
    sc.pp.filter_genes(ad,min_cells=3)
    if ad.n_obs<40 or ad.n_vars<100: raise RuntimeError('Insufficient cells/genes after processing')
    sc.pp.normalize_total(ad,target_sum=1e4)
    sc.pp.log1p(ad)
    ad.layers['lognorm']=ad.X.copy()
    n_top=min(2500,max(100,ad.n_vars-1))

    sc.pp.highly_variable_genes(ad,n_top_genes=n_top,flavor='seurat',batch_key='sample_id')
    n_pcs=min(args.n_pcs,int(ad.var['highly_variable'].sum())-1,ad.n_obs-2)
    if n_pcs<5: raise RuntimeError('Too few variable features for PCA')
    sc.pp.pca(ad,n_comps=n_pcs,mask_var='highly_variable',svd_solver='arpack',random_state=1234)
    sc.pp.neighbors(ad,n_neighbors=15,n_pcs=n_pcs,random_state=1234)
    sc.tl.umap(ad,random_state=1234)
    sc.tl.leiden(ad,resolution=args.resolution,flavor='igraph',directed=False,random_state=1234)
    status(registry,'normalization_clustering','COMPLETE',
           f'Library-size normalized log1p; HVGs; PCA {n_pcs}; Leiden {args.resolution}')
    status(registry,'integration','NOT_APPLIED',
           'No forced batch correction. Inspect clusters that are dominated by one donor before comparing conditions.')
    for name,label in [('leiden','clusters'),('sample_id','samples'),('cell_type','putative_celltypes')]:
        if name=='cell_type':continue
        umap_plot(ad,name,figures,f'UMAP_{label}')
    clonotypes(ad,inp,samples,tables,registry)
    annotate(ad,args.markers,tables,registry)
    tcr_guard(ad,tables,registry)
    cluster_markers(ad,tables,figures,registry)
    umap_plot(ad,'cell_type',figures,'UMAP_putative_celltypes')
    immune_subclustering(ad,tables,figures,registry)
    for gene in ['CD3D','CD8A','MS4A1','LYZ','EPCAM','KRT19','C1QA','PDCD1','GZMB','FOXP3']:
        feature_umap(ad,gene,figures)
    comp=pd.crosstab(ad.obs['sample_id'],ad.obs['cell_type'])
    compositions(ad,tables,registry)
    found=signature_scores(ad,analysis_cfg['scoring']['immune_modules'],
                           int(analysis_cfg['scoring'].get('min_genes_present',2)),'signature')
    pathway_plot(ad,list(found),figures)
    for module in ['Cytotoxicity','Exhaustion_associated','Treg_associated','Myeloid_SPP1','IFNG_response']:
        signature_umap(ad,module,figures)
    status(registry,'immune_states','EXPLORATORY',
           f'{len(found)} RNA expression programs. No functional exhaustion/protein proof')
    pseudobulk(ad,tables,registry)
    ligand_receptor(ad,args.ligand_receptor,tables,registry)
    trajectory(ad,tables,registry)
    targets(ad,tables,registry)
    status(registry,'malignant_CNV_classification','NOT_RUN',
           'Require valid gene coordinate reference and appropriate normal baseline; epithelial markers are not malignancy proof')
    status(registry,'ambient_RNA_correction','NOT_RUN',
           'Unfiltered background assessed; correction requires validated SoupX/CellBender model')
    status(registry,'reference_annotation','NOT_RUN',
           'Marker-based preliminary annotations only; external cancer-specific references required')
    status(registry,'TF_regulons_SCENIC','NOT_RUN',
           'Transcription factor expression is not TF activity; validated regulon reference required')
    status(registry,'NicheNet_LIANA_CellChat','NOT_RUN',
           'Only candidate ligand/receptor co-expression module executed; full inference not executed')
    status(registry,'BCR_CITEseq','NOT_AVAILABLE',
           'No BCR or ADT files were provided')
    status(registry,'spatial_immune_exclusion','NOT_AVAILABLE',
           'No spatial coordinates with dissociated tumor scRNA-seq')
    status(registry,'response_ML_survival','NOT_TESTABLE',
           'No responder labels or survival outcomes are in the sample sheet')
    status(registry,'external_validation','NOT_RUN','Needs independent cohort and clinical hypothesis')
    status(registry,'drug_discovery','NOT_RUN',
           'Target expression hypotheses only, no structure screening or functional evidence')
    status(registry,'experimental_validation','NOT_AVAILABLE','Requires wet lab assays')
    import atac
    atac.run(ad, samples, inp, out, figures, tables, registry, analysis_cfg,
             args.gene_annotation or None)
    provenance={'pipeline':'scRNA-seq-scATAC-seq-end2end','input_samples':samples.sample_id.tolist(),
                'source':'Sample sheet passed with --samples. For this cohort the tumor-type column is the GEO characteristic tumor type.',
                'scanpy_version':sc.__version__,'python':sys.version,
                'parameters':vars(args),'seed':1234,
                'caution':'Exploratory single-cell atlas, not publication-validated claims'}
    (out/'provenance.json').write_text(json.dumps(provenance,indent=2))
    overview=pd.DataFrame(overview)
    report(out,registry,samples,overview,ad)
    LOG.info('SUCCESS: %d cells %d genes, %d figure PDFs',ad.n_obs,ad.n_vars,len(list(figures.glob('*.pdf'))))


if __name__=='__main__':main()
