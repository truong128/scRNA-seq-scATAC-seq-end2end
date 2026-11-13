# Pipeline

The repository README at the top level is the guide for this end-to-end single-cell RNA-seq and single-cell ATAC-seq workflow. It covers Ubuntu and macOS.

From the repository root:

```bash
./run_end2end.sh
```

`main.nf` is the Nextflow workflow. `bin/analyze.py` runs single-cell RNA-seq, TCR, and the call into single-cell ATAC-seq. `bin/atac.py` runs ATAC QC, LSI, and gene activity. `bin/edger_pseudobulk.R` runs the donor-level edgeR test. `bin/verify_local_inputs.py` checks `inputs/GEO_SOURCE.tsv` before an offline rerun.
