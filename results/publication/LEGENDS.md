# Legends

Citation: Truong Nguyen, CUNY School of Medicine.

Data: GEO series GSE274934. Single-cell RNA-seq and single-cell ATAC-seq are separate libraries. Cells are not matched across the two assays.

## Figure 1

`figures/Figure1_single_cell_rna_and_atac.png`

Single-cell RNA-seq and single-cell ATAC-seq from the same donors, analysed as separate libraries.

**A.** UMAP of 97,017 single-cell RNA-seq cells that passed QC (200 to 9,000 genes, mitochondrial fraction at most 25%). Colors are marker-based labels. `NK_TCRpos_ambiguous`, `B_cells_TCRpos_ambiguous`, and `Unresolved` are shown because those cells remain in the object. They are not confirmed cell types, and they are excluded from ligand–receptor pairs, the target table, and edgeR. No batch correction was applied. The panel is exploratory.

**B.** UMAP of gene activity for 24,712 single-cell ATAC-seq cells from AL05, AL08, AL10, AL11, and AL12. Colors are donor libraries. Gene activity counts GENCODE v44 basic protein-coding and lncRNA genes, longest transcript, TSS ± 2,000 bp. Clusters in this panel are not RNA cell types. Peak sets were not concatenated, and no same-cell model was run.

Vector copies: `figures/Figure1A_scrna_umap.pdf` and `figures/Figure1B_scatac_gene_activity_umap.pdf`.

Supporting PDF and PNG files, including QC, markers, per-library ATAC LSI, and the donor-level RNA versus gene-activity comparison, are in `results/analysis/figures/`.

## Tables

**Table 1.** `tables/Table1_cohort.csv`. Donors, tumor type from the GEO characteristic `tumor type`, which single-cell assays are present, and GEO sample accessions. A blank ATAC accession means that donor has no ATAC library.

**Table 2.** `tables/Table2_scrna_qc.csv`. Single-cell RNA-seq input cells, retained cells, and median genes and UMIs per donor.

**Table 3.** `tables/Table3_scatac_qc.csv`. Single-cell ATAC-seq cells before and after the pre-set QC, fragment summaries, and whether LSI component 1 was dropped because it tracked depth.

**Table 4.** `tables/Table4_scrna_scatac_barcode_overlap.csv`. Shared barcodes between the single-cell RNA-seq library and the single-cell ATAC-seq library of the same donor. The largest shared fraction is 0.0004.

**Table 5.** `tables/Table5_scrna_scatac_sample_correlation.csv`. Donor-level Spearman correlation of mean RNA expression and mean ATAC gene activity. The comparison is between separate libraries.

**Table 6.** `tables/Table6_edger_summary.csv`. Unpaired donor-level edgeR quasi-likelihood F test. Contrast: `WT_LUAD vs ALK_positive_LUAD`. A positive logFC means higher in WT. `FDR_0_05 = 0` means the test ran and no gene passed FDR 0.05.

**Table 7.** `tables/Table7_labels_not_tested.csv`. Labels excluded from edgeR. `Unresolved` is not a cell type. A `TCRpos_ambiguous` label is not a confirmed cell type.

**Table 8.** `tables/Table8_cell_counts.csv`. Cells per donor and marker label, including the labels that were not tested.
