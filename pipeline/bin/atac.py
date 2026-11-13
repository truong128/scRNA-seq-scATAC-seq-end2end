#!/usr/bin/env python3

from __future__ import annotations

import gzip
from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
from scipy import sparse, stats
from sklearn.decomposition import TruncatedSVD
import scanpy as sc

import analyze as A


def _cfg(analysis_cfg):
    raw = (analysis_cfg or {}).get('atac') or {}
    return {
        'min_fragments': int(raw.get('min_fragments', 1000)),
        'min_fraction_in_peaks': float(raw.get('min_fraction_in_peaks', 0.20)),
        'min_fraction_tss_fragments': float(raw.get('min_fraction_tss_fragments', 0.15)),
        'promoter_window_bp': int(raw.get('promoter_window_bp', 2000)),
        'lsi_components': int(raw.get('lsi_components', 30)),
        'drop_lsi1_if_depth_correlation': float(raw.get('drop_lsi1_if_depth_correlation', 0.5)),
        'leiden_resolution': float(raw.get('leiden_resolution', 0.6)),
        'gene_types': tuple(raw.get('gene_types', ['protein_coding', 'lncRNA'])),
    }


def find_atac_samples(folder: Path, sample_ids):
    found = []
    for sid in sample_ids:
        peaks = folder / f'{sid}_atac_peaks.h5'
        single = folder / f'{sid}_atac_singlecell.csv.gz'
        if peaks.exists() and single.exists():
            found.append(sid)
        elif peaks.exists() or single.exists():
            raise FileNotFoundError(f'{sid} has only one of the two ATAC files ({peaks.name}, {single.name})')
    return found


def load_promoters(gtf: Path, window: int, gene_types):
    best = {}
    opener = gzip.open if str(gtf).endswith('.gz') else open
    with opener(gtf, 'rt') as handle:
        for line in handle:
            if not line or line[0] == '#':
                continue
            parts = line.rstrip('\n').split('\t')
            if len(parts) < 9 or parts[2] != 'transcript':
                continue
            chrom, start_s, end_s, strand, attrs = parts[0], parts[3], parts[4], parts[6], parts[8]
            if not chrom.startswith('chr') or chrom in ('chrM', 'chrMT'):
                continue
            info = {}
            for field in attrs.strip().split(';'):
                field = field.strip()
                if not field or ' ' not in field:
                    continue
                key, value = field.split(' ', 1)
                info[key] = value.strip().strip('"')
            if info.get('gene_type') not in gene_types:
                continue
            gene = info.get('gene_name') or info.get('gene_id')
            start, end = int(start_s), int(end_s)
            length = end - start + 1
            previous = best.get(gene)
            if previous is None or length > previous[0]:
                tss = start if strand == '+' else end
                best[gene] = (length, chrom, tss)
    rows = []
    for gene, (_, chrom, tss) in best.items():

        rows.append((chrom, tss - 1 - window, tss + window, gene))
    out = pd.DataFrame(rows, columns=['chrom', 'start0', 'end0', 'gene'])
    if out.empty:
        raise RuntimeError(f'No promoters parsed from {gtf}')
    return out


def load_peaks(path: Path):
    import h5py
    with h5py.File(path, 'r') as handle:
        shape = handle['matrix/shape'][:]
        matrix = sparse.csc_matrix(
            (handle['matrix/data'][:], handle['matrix/indices'][:], handle['matrix/indptr'][:]),
            shape=tuple(int(x) for x in shape),
        ).T.tocsr()
        barcodes = [x.decode() for x in handle['matrix/barcodes'][:]]
        names = [x.decode() for x in handle['matrix/features/name'][:]]
        kind = {x.decode() for x in handle['matrix/features/feature_type'][:]}
    if kind != {'Peaks'}:
        raise RuntimeError(f'{path.name} feature types are {kind}, expected only Peaks')
    chrom, start, end = [], [], []
    for name in names:
        left, right = name.split(':')
        a, b = right.split('-')
        chrom.append(left)
        start.append(int(a))
        end.append(int(b))
    peaks = pd.DataFrame({'chrom': chrom, 'start0': start, 'end0': end})
    return matrix, barcodes, peaks


def cell_qc(singlecell: Path, barcodes, cfg):
    wanted = set(barcodes)
    rows = []
    with gzip.open(singlecell, 'rt') as handle:
        for row in pd.read_csv(handle, chunksize=200_000):
            hit = row[row['barcode'].isin(wanted)]
            if len(hit):
                rows.append(hit)
    if not rows:
        raise RuntimeError(f'None of the matrix barcodes are in {singlecell.name}')
    table = pd.concat(rows, ignore_index=True).drop_duplicates('barcode')
    table = table.set_index('barcode').reindex(barcodes)
    if table['passed_filters'].isna().any():
        missing = int(table['passed_filters'].isna().sum())
        raise RuntimeError(f'{missing} matrix barcodes have no singlecell QC row')
    passed = table['passed_filters'].astype(float)
    table['fraction_tss_fragments'] = table['TSS_fragments'].astype(float) / passed.replace(0, np.nan)
    table['fraction_in_peaks'] = table['peak_region_fragments'].astype(float) / passed.replace(0, np.nan)
    keep = (
        (passed >= cfg['min_fragments'])
        & (table['fraction_in_peaks'] >= cfg['min_fraction_in_peaks'])
        & (table['fraction_tss_fragments'] >= cfg['min_fraction_tss_fragments'])
    ).fillna(False).to_numpy()
    return table, keep


def gene_activity(matrix, peaks, promoters):
    gene_names = promoters['gene'].drop_duplicates().tolist()
    gene_index = {gene: i for i, gene in enumerate(gene_names)}
    peak_i, gene_i = [], []
    for chrom, peak_sub in peaks.groupby('chrom', sort=False):
        window = promoters[promoters['chrom'] == chrom]
        if window.empty or peak_sub.empty:
            continue
        starts = peak_sub['start0'].to_numpy()
        ends = peak_sub['end0'].to_numpy()
        indices = peak_sub.index.to_numpy()
        overlap = (starts[:, None] < window['end0'].to_numpy()[None, :]) & (
            ends[:, None] > window['start0'].to_numpy()[None, :]
        )
        rr, cc = np.nonzero(overlap)
        peak_i.append(indices[rr])
        gene_i.append(np.array([gene_index[g] for g in window['gene'].to_numpy()[cc]], dtype=np.int32))
    if not peak_i:
        raise RuntimeError('No peak overlapped a promoter window')
    mapping = sparse.csr_matrix(
        (np.ones(sum(len(x) for x in peak_i), dtype=np.float32),
         (np.concatenate(peak_i), np.concatenate(gene_i))),
        shape=(matrix.shape[1], len(gene_names)),
    )
    activity = matrix @ mapping
    return activity.tocsr(), gene_names, int(mapping.nnz)


def lsi(binary, n_components, depth, drop_at):
    n_components = int(min(n_components, binary.shape[0] - 1, binary.shape[1] - 1))
    if n_components < 2:
        raise RuntimeError(f'Only {binary.shape[0]} cells / {binary.shape[1]} peaks; LSI needs more')
    present = np.asarray((binary > 0).sum(axis=0)).ravel()
    keep = present >= max(5, int(0.01 * binary.shape[0]))
    if keep.sum() < 50:
        raise RuntimeError('Fewer than 50 peaks pass the 1% cell filter')
    used = binary[:, keep]
    idf = np.log1p(used.shape[0] / (1.0 + present[keep]))
    weighted = used.multiply(idf).tocsr()
    depths = np.asarray(weighted.sum(axis=1)).ravel()
    depths[depths == 0] = 1
    weighted = sparse.diags(1.0 / depths) @ weighted
    svd = TruncatedSVD(n_components=n_components, random_state=1234)
    scores = svd.fit_transform(weighted)
    corr = float(stats.spearmanr(scores[:, 0], np.log10(depth + 1)).correlation)
    drop = abs(corr) >= drop_at
    used_scores = scores[:, 1:] if drop and scores.shape[1] > 2 else scores
    return used_scores, {'lsi1_depth_spearman': round(corr, 4), 'lsi1_dropped': bool(drop),
                         'peaks_used': int(keep.sum()), 'components_used': int(used_scores.shape[1])}


def _umap(ad, color, path_stem, title, note):
    fig, ax = plt.subplots(figsize=(7.2, 5.4))
    xy = ad.obsm['X_umap']
    if ad.obs[color].dtype.name == 'category' or ad.obs[color].dtype == object:
        labels = ad.obs[color].astype(str)
        for lab in pd.unique(labels):
            mask = (labels == lab).to_numpy()
            ax.scatter(xy[mask, 0], xy[mask, 1], s=3, linewidths=0, label=lab, rasterized=True)
        ax.legend(markerscale=3, frameon=False, fontsize=7, loc='center left', bbox_to_anchor=(1.02, 0.5))
    else:
        points = ax.scatter(xy[:, 0], xy[:, 1], c=ad.obs[color], s=3, linewidths=0, rasterized=True, cmap='viridis')
        fig.colorbar(points, ax=ax, label=color)
    ax.set(title=title, xlabel='UMAP 1', ylabel='UMAP 2')
    A.stamp(ax, note)
    A.savefig(fig, path_stem.parent, path_stem.name)


def run(rna, samples, matrix_dir, outdir, figures, tables, registry, analysis_cfg, gtf):
    cfg = _cfg(analysis_cfg)
    matrix_dir = Path(matrix_dir)
    sample_ids = [str(x) for x in samples['sample_id'].tolist()]
    atac_ids = find_atac_samples(matrix_dir, sample_ids)
    if not atac_ids:
        A.status(registry, 'ATAC', 'NOT_AVAILABLE',
                 'No SAMPLE_atac_peaks.h5 + SAMPLE_atac_singlecell.csv.gz next to the RNA matrices')
        return
    if gtf is None or not Path(gtf).exists():
        A.status(registry, 'ATAC_gene_activity', 'FAILED',
                 'GENCODE GTF missing; per-sample QC will still run but gene activity and RNA linking will not')
        gtf = None
    else:
        promoters = load_promoters(Path(gtf), cfg['promoter_window_bp'], cfg['gene_types'])
        A.status(registry, 'ATAC_gene_annotation', 'COMPLETE',
                 f'{promoters.gene.nunique()} genes, longest-transcript TSS, '
                 f'+/- {cfg["promoter_window_bp"]} bp, GENCODE file {Path(gtf).name}, types {cfg["gene_types"]}')

    meta = samples.set_index('sample_id')
    activities = []
    qc_rows = []
    overlap_rows = []
    for sid in atac_ids:
        matrix, barcodes, peaks = load_peaks(matrix_dir / f'{sid}_atac_peaks.h5')
        table, keep = cell_qc(matrix_dir / f'{sid}_atac_singlecell.csv.gz', barcodes, cfg)
        rna_barcodes = set(rna.obs.loc[rna.obs['sample_id'].astype(str) == sid, 'orig_barcode'].astype(str))
        shared = len(rna_barcodes & set(barcodes))
        overlap_rows.append({
            'sample_id': sid, 'rna_qc_cells': len(rna_barcodes), 'atac_matrix_cells': len(barcodes),
            'shared_barcodes': shared,
            'shared_fraction_of_smaller': round(shared / max(1, min(len(rna_barcodes), len(barcodes))), 6),
        })
        n_in, n_keep = len(barcodes), int(keep.sum())
        A.status(registry, f'{sid}:ATAC_QC', 'COMPLETE',
                 f'{n_keep}/{n_in} cells kept. Filters: fragments>={cfg["min_fragments"]}, '
                 f'fraction in peaks>={cfg["min_fraction_in_peaks"]}, '
                 f'fraction of fragments overlapping the Cell Ranger TSS region>={cfg["min_fraction_tss_fragments"]}. '
                 'This fraction is not a TSS enrichment score.')
        if n_keep < 50:
            A.status(registry, f'{sid}:ATAC_LSI', 'SKIPPED', f'Only {n_keep} cells after QC')
            continue
        kept = matrix[keep]
        depth = table['passed_filters'].to_numpy()[keep]
        binary = kept.copy()
        binary.data = np.ones_like(binary.data)
        try:
            scores, info = lsi(binary, cfg['lsi_components'], depth, cfg['drop_lsi1_if_depth_correlation'])
        except RuntimeError as exc:
            A.status(registry, f'{sid}:ATAC_LSI', 'FAILED', str(exc))
            continue
        obs = pd.DataFrame(index=[f'{sid}_{barcodes[i]}' for i, flag in enumerate(keep) if flag])
        obs['sample_id'] = sid
        obs['patient_id'] = str(meta.loc[sid, 'patient_id'])
        obs['condition'] = str(meta.loc[sid, 'condition'])
        obs['passed_filters'] = depth
        obs['fraction_in_peaks'] = table['fraction_in_peaks'].to_numpy()[keep]
        obs['fraction_tss_fragments'] = table['fraction_tss_fragments'].to_numpy()[keep]

        one = sc.AnnData(sparse.csr_matrix((len(obs), 1), dtype=np.float32), obs=obs)
        one.obsm['X_lsi'] = scores
        sc.pp.neighbors(one, use_rep='X_lsi', n_neighbors=15, random_state=1234)
        sc.tl.leiden(one, resolution=cfg['leiden_resolution'], flavor='igraph', directed=False, random_state=1234)
        sc.tl.umap(one, random_state=1234)
        _umap(one, 'leiden', figures / f'ATAC_LSI_{sid}',
              f'{sid} ATAC LSI (own peaks)',
              'Within-sample peak LSI. Not comparable across samples: each library has its own peak set.')
        A.status(registry, f'{sid}:ATAC_LSI', 'COMPLETE',
                 f'Binarized TF-IDF LSI. {info}. Leiden {cfg["leiden_resolution"]}. '
                 f'LSI1 dropped when |Spearman| with log10(fragments) >= {cfg["drop_lsi1_if_depth_correlation"]}.')
        one.obs[['sample_id', 'leiden', 'passed_filters', 'fraction_in_peaks', 'fraction_tss_fragments']].to_csv(
            tables / f'ATAC_{sid}_lsi_cells.csv')
        qc_rows.append({
            'sample_id': sid, 'condition': obs['condition'].iloc[0], 'cells_in_matrix': n_in,
            'cells_after_qc': n_keep, 'median_fragments': float(np.median(depth)),
            'median_fraction_in_peaks': float(np.median(obs['fraction_in_peaks'])),
            'median_fraction_tss_fragments': float(np.median(obs['fraction_tss_fragments'])),
            **info,
        })
        if gtf is not None:
            activity, genes, n_links = gene_activity(kept, peaks, promoters)
            act = sc.AnnData(activity, obs=obs.copy())
            act.var_names = genes
            activities.append(act)
            A.status(registry, f'{sid}:ATAC_gene_activity', 'COMPLETE',
                     f'{activity.shape[1]} genes, {n_links} peak-gene window overlaps. '
                     'A peak is counted for every gene it overlaps.')
        del matrix, kept, binary, one

    pd.DataFrame(overlap_rows).to_csv(tables / 'ATAC_RNA_barcode_overlap.csv', index=False)
    pd.DataFrame(qc_rows).to_csv(tables / 'ATAC_sample_qc.csv', index=False)
    max_shared = max((row['shared_fraction_of_smaller'] for row in overlap_rows), default=0)
    if max_shared < 0.01:
        A.status(registry, 'ATAC_RNA_cell_integration', 'NOT_APPLICABLE',
                 f'Separate libraries. The largest shared-barcode fraction is {max_shared:.4f}. '
                 'No WNN, MultiVI, or same-cell peak-to-gene link was run. See ATAC_RNA_barcode_overlap.csv.')
    else:
        A.status(registry, 'ATAC_RNA_cell_integration', 'NOT_RUN',
                 f'Shared barcodes reach {max_shared:.4f}. Cell-level integration is not implemented; '
                 'do not treat the sample-level correlation below as a paired-cell result.')

    donors = pd.DataFrame(qc_rows)
    if len(donors):
        counts = donors.groupby('condition')['sample_id'].nunique()
        if len(counts) < 2 or (counts < 3).any():
            A.status(registry, 'ATAC_between_condition', 'NOT_TESTABLE',
                     f'No differential accessibility test. ATAC donors by tumor type: {counts.to_dict()}. '
                     'Need >=3 donors in every group. RNA can still be tested if RNA replication is sufficient.')
        else:
            A.status(registry, 'ATAC_between_condition', 'NOT_RUN',
                     'Replication would allow a donor-level test, but no ATAC differential model is implemented.')

    if gtf is None or not activities:
        return
    joint = sc.concat(activities, join='inner', index_unique=None)
    del activities
    sc.pp.filter_genes(joint, min_cells=20)
    sc.pp.normalize_total(joint, target_sum=1e4)
    sc.pp.log1p(joint)
    n_top = min(2000, joint.n_vars - 1)
    sc.pp.highly_variable_genes(joint, n_top_genes=n_top, flavor='seurat', batch_key='sample_id')
    n_pcs = min(30, int(joint.var['highly_variable'].sum()) - 1, joint.n_obs - 2)
    sc.pp.pca(joint, n_comps=n_pcs, mask_var='highly_variable', svd_solver='arpack', random_state=1234)
    sc.pp.neighbors(joint, n_neighbors=15, n_pcs=n_pcs, random_state=1234)
    sc.tl.leiden(joint, resolution=cfg['leiden_resolution'], flavor='igraph', directed=False, random_state=1234)
    sc.tl.umap(joint, random_state=1234)
    note = 'Gene-activity embedding across separate ATAC libraries. Clusters are not RNA cell types.'
    _umap(joint, 'sample_id', figures / 'ATAC_gene_activity_UMAP_samples', 'ATAC gene activity by sample', note)
    _umap(joint, 'leiden', figures / 'ATAC_gene_activity_UMAP_clusters', 'ATAC gene activity clusters', note)
    joint.obs[['sample_id', 'patient_id', 'condition', 'leiden']].to_csv(tables / 'ATAC_gene_activity_cells.csv')
    purity = joint.obs.groupby('leiden')['sample_id'].value_counts(normalize=True).groupby(level=0).max()
    A.status(registry, 'ATAC_joint_embedding', 'EXPLORATORY',
             f'{joint.n_obs} cells, {joint.n_vars} genes after the min-20-cells filter. '
             f'{int((purity > 0.9).sum())}/{len(purity)} clusters are >90% one sample. No batch correction was applied.')

    rows = []
    for sid in joint.obs['sample_id'].astype(str).unique():
        rna_sub = rna[rna.obs['sample_id'].astype(str) == sid]
        atac_sub = joint[joint.obs['sample_id'].astype(str) == sid]
        if rna_sub.n_obs < 30 or atac_sub.n_obs < 30:
            continue
        genes = sorted(set(rna_sub.var_names) & set(atac_sub.var_names))
        rna_mean = np.asarray(rna_sub[:, genes].X.mean(axis=0)).ravel()
        atac_mean = np.asarray(atac_sub[:, genes].X.mean(axis=0)).ravel()
        use = (rna_mean > 0) & (atac_mean > 0)
        if use.sum() < 200:
            rows.append({'sample_id': sid, 'n_genes': int(use.sum()), 'spearman': np.nan})
            continue
        rho = float(stats.spearmanr(rna_mean[use], atac_mean[use]).correlation)
        rows.append({'sample_id': sid, 'condition': str(meta.loc[sid, 'condition']),
                     'n_genes': int(use.sum()), 'spearman_rna_vs_gene_activity': round(rho, 4),
                     'rna_cells': int(rna_sub.n_obs), 'atac_cells': int(atac_sub.n_obs)})
    corr = pd.DataFrame(rows)
    corr.to_csv(tables / 'ATAC_RNA_sample_correlation.csv', index=False)
    if len(corr) and corr.get('spearman_rna_vs_gene_activity') is not None and corr['spearman_rna_vs_gene_activity'].notna().any():
        fig, ax = plt.subplots(figsize=(6.5, 4))
        ax.bar(corr['sample_id'], corr['spearman_rna_vs_gene_activity'], color='#5D8FA6')
        ax.set(ylabel='Spearman', ylim=(0, 1), title='Donor-level RNA vs ATAC gene activity')
        A.stamp(ax, 'Same donor, different cells. Confounded by cell composition. Not a paired-cell correlation.')
        A.savefig(fig, figures, 'ATAC_RNA_sample_correlation')
        A.status(registry, 'ATAC_RNA_sample_integration', 'SAMPLE_LEVEL_ONLY',
                 'Spearman of mean RNA expression vs mean gene activity within each donor that has both assays. '
                 'Cells are not paired. Composition differs between the RNA and ATAC libraries. '
                 f'Values: ' + ', '.join(f"{r.sample_id}={r.spearman_rna_vs_gene_activity}" for r in corr.itertuples()))
    else:
        A.status(registry, 'ATAC_RNA_sample_integration', 'FAILED', 'No donor had 200 genes measured in both assays')
