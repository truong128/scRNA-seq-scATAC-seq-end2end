#!/usr/bin/env python3

import argparse
import gzip
from pathlib import Path
import numpy as np
from scipy import sparse
import h5py

GENES = ['CD3D','CD3E','TRAC','CD2','CD4','IL7R','CD8A','CD8B','GZMK','FOXP3',
         'IL2RA','CTLA4','TIGIT','PDCD1','LAG3','HAVCR2','TOX','ENTPD1','CXCL13',
         'NKG7','GNLY','KLRD1','FCGR3A','NCR1','MS4A1','CD79A','CD79B','CD37',
         'JCHAIN','MZB1','XBP1','SDC1','LYZ','FCN1','S100A8','S100A9','C1QA',
         'C1QB','C1QC','APOE','CD68','FCER1A','CD1C','CLEC10A','HLA-DRA',
         'CLEC4C','IL3RA','LILRA4','KIT','TPSAB1','TPSB2','COL1A1','COL1A2',
         'DCN','LUM','PECAM1','VWF','KDR','CLDN5','EPCAM','KRT8','KRT18',
         'KRT19','PRF1','GZMB','GZMA','STAT1','IRF1','GBP1','CXCL9','CXCL10',
         'IDO1','HLA-A','HLA-B','HLA-C','B2M','TAP1','TAP2','NLRC5','SPP1',
         'TREM2','LPL','VEGFA','CA9','SLC2A1','LDHA','PDK1','MKI67','TOP2A',
         'STMN1','PCNA','FOS','JUN','JUNB','ATF3','DUSP1','HSPA1A','CD274',
         'PDCD1LG2','PVR','CD226','CD80','CD86','CD28','CXCR3','CXCR4','CCL2',
         'CCR2','CSF1','CSF1R','IL10','IL10RA','TGFB1','TGFBR2','IFNG',
         'IFNGR1','TNF','TNFRSF1A','CD44','MIF','CD74','LGALS9','JAK1',
         'JAK2','STAT3','CD47','TCF7','SLAMF6','HSPA1B','MT-CO1','MT-ND1']
GENES = list(dict.fromkeys(GENES)) + [f'SYNTHETIC_GENE_{i:04d}' for i in range(700)]


def put_10x(path, X, barcodes, genes):

    csc=sparse.csc_matrix(X.T,dtype=np.int32)
    with h5py.File(path,'w') as f:
        g=f.create_group('matrix')
        for key,v in [('data',csc.data),('indices',csc.indices),
                      ('indptr',csc.indptr),('shape',np.asarray(csc.shape,dtype=np.int64))]:
            g.create_dataset(key,data=v)
        g.create_dataset('barcodes',data=np.asarray(barcodes,dtype='S'))
        feat=g.create_group('features')
        feat.create_dataset('id',data=np.asarray(['syn_'+g for g in genes],dtype='S'))
        feat.create_dataset('name',data=np.asarray(genes,dtype='S'))
        feat.create_dataset('feature_type',data=np.asarray(['Gene Expression']*len(genes),dtype='S'))
        feat.create_dataset('genome',data=np.asarray(['GRCh38']*len(genes),dtype='S'))
        feat.create_dataset('_all_tag_keys',data=np.asarray(['genome'],dtype='S'))


def main():
    p=argparse.ArgumentParser()
    p.add_argument('--out',default='tests/synthetic_matrices')
    args=p.parse_args()
    out=Path(args.out)
    out.mkdir(parents=True,exist_ok=True)
    rng=np.random.default_rng(42)
    n=220
    for sample in ['TEST01','TEST02']:
        counts=rng.poisson(.30,(n,len(GENES))).astype(np.int16)
        lineages=['T','Macrophages','Epithelial','B_cells']
        programs={
            'T':['CD3D','CD3E','TRAC','CD2','CD8A','CD8B'],
            'Macrophages':['LYZ','C1QA','C1QB','C1QC','CD68','APOE'],
            'Epithelial':['EPCAM','KRT8','KRT18','KRT19'],
            'B_cells':['MS4A1','CD79A','CD79B','CD37'],
        }
        for i in range(n):
            typ=lineages[i%len(lineages)]
            for name in programs[typ]:counts[i,GENES.index(name)]+=rng.poisson(8)
        barcodes=[f'AAAACCCCGGGG{i:04d}-1' for i in range(n)]
        put_10x(out/f'{sample}_filtered.h5',counts,barcodes,GENES)
        background=rng.poisson(.002,(100,len(GENES))).astype(np.int16)
        put_10x(out/f'{sample}_raw.h5',np.vstack([counts,background]),
                barcodes+[f'BBBBCCCCGGGG{i:04d}-1' for i in range(100)],GENES)
        with gzip.open(out/f'{sample}_tcr.csv.gz','wt') as f:
            f.write('barcode,chain,productive,raw_clonotype_id\n')
            for i in range(0,n,4):
                f.write(f'{barcodes[i]},TRA,True,{sample}_clone{i%20}\n')
    (out/'synthetic_samples.csv').write_text(
        'sample_id,patient_id,condition\nTEST01,TEST01,SYNTHETIC\nTEST02,TEST02,SYNTHETIC\n')
    print('Synthetic fixture created in',out,'DO NOT interpret as real biological data.')

if __name__=='__main__':main()
