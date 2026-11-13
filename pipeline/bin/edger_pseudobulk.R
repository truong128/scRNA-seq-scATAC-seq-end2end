suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly=TRUE)
if (length(args)<2) stop('Usage: Rscript edger_pseudobulk.R analysis/tables output_folder')
in_dir <- args[[1]]
out_dir <- args[[2]]
dir.create(out_dir, showWarnings=FALSE, recursive=TRUE)
meta <- read.csv(file.path(in_dir,'pseudobulk_metadata.csv'), stringsAsFactors=FALSE)
if (!nrow(meta)) stop('No pseudobulk meta')
cond <- unique(meta[,c('patient_id','condition')])
replicates <- table(cond$condition)
if (length(replicates)<2 || any(replicates<3)) {
  writeLines(sprintf('NOT TESTABLE: at least 3 independent donors per condition recommended; observed %s',
    paste(names(replicates),as.integer(replicates),collapse='; ')),file.path(out_dir,'NOT_TESTABLE.txt'))
  quit(status=0)
}
if (anyDuplicated(cond$patient_id)) {
  stop('Repeated patient with multiple conditions: use paired design, do not use this unpaired script')
}
mat <- as.matrix(read.csv(gzfile(file.path(in_dir,'pseudobulk_gene_by_sample_celltype.csv.gz')),
                          row.names=1,check.names=FALSE))
stopifnot(all(meta$pseudobulk_id %in% colnames(mat)))
if (length(unique(meta$condition))!=2) stop('Only exactly two conditions in this conservative template')
levels <- sort(unique(meta$condition))
contrast <- paste(levels[2],'vs',levels[1])
results_summary <- list()
skipped <- list()
for (ctype in unique(meta$cell_type)) {
  smeta <- meta[meta$cell_type==ctype,,drop=FALSE]
  donor_cond <- unique(smeta[,c('patient_id','condition')])
  reps <- table(donor_cond$condition)
  reps_full <- setNames(integer(length(levels)), levels)
  reps_full[names(reps)] <- as.integer(reps)
  donor_txt <- paste(names(reps_full), reps_full, sep='=', collapse='; ')
  if (ctype == 'Unresolved' || grepl('TCRpos_ambiguous', ctype)) {
    reason <- if (ctype == 'Unresolved') {
      'Unresolved is not a cell type, so it is not tested'
    } else {
      'TCR-ambiguous label is not a confirmed cell type, so it is not tested'
    }
    skipped[[length(skipped)+1]] <- data.frame(
      celltype=ctype, reason=reason, donors=donor_txt, stringsAsFactors=FALSE)
    next
  }
  if (length(reps)!=2 || any(reps<3)) {
    skipped[[length(skipped)+1]] <- data.frame(
      celltype=ctype, reason='fewer than 3 independent donors in each condition',
      donors=donor_txt, stringsAsFactors=FALSE)
    next
  }
  if (nrow(smeta) != nrow(donor_cond)) {
    stop(sprintf(paste0('Cell type %s has several samples per donor. They are not independent replicates: ',
      'aggregate per donor or use duplicateCorrelation/a mixed model before running this script.'), ctype))
  }
  cts <- round(mat[,smeta$pseudobulk_id,drop=FALSE])
  colnames(cts) <- smeta$pseudobulk_id
  if (nrow(cts)<10 || ncol(cts)<6) {
    skipped[[length(skipped)+1]] <- data.frame(
      celltype=ctype, reason='fewer than 6 donor libraries or 10 genes',
      donors=donor_txt, stringsAsFactors=FALSE)
    next
  }
  design <- model.matrix(~ condition,data=transform(smeta,condition=factor(condition,levels=levels)))
  if (qr(design)$rank<ncol(design)) stop('Design not full rank')
  y <- DGEList(cts)
  keep <- filterByExpr(y,design=design)
  if (sum(keep)<100) {
    skipped[[length(skipped)+1]] <- data.frame(
      celltype=ctype, reason='filterByExpr kept fewer than 100 genes',
      donors=donor_txt, stringsAsFactors=FALSE)
    next
  }
  y <- y[keep,,keep.lib.sizes=FALSE]
  y <- calcNormFactors(y)
  y <- estimateDisp(y,design,robust=TRUE)
  fit <- glmQLFit(y,design,robust=TRUE)
  res <- topTags(glmQLFTest(fit,coef=2),n=Inf)$table
  res$gene <- rownames(res)
  res$contrast <- contrast
  res$celltype <- ctype
  target <- gsub('[^A-Za-z0-9_]+','_',ctype)
  write.csv(res,file.path(out_dir,paste0('DE_',target,'.csv')),row.names=FALSE)
  results_summary[[ctype]] <- data.frame(celltype=ctype,tested_genes=nrow(res),
                                         FDR_0_05=sum(res$FDR<.05),contrast=contrast,
                                         stringsAsFactors=FALSE)
}
if (length(results_summary)) do.call(rbind,results_summary) |>
  write.csv(file.path(out_dir,'DE_summary.csv'),row.names=FALSE)
if (!length(results_summary)) writeLines('No cell types with >=3 donors per condition and sufficient counts',
                                         file.path(out_dir,'NOT_TESTABLE.txt'))
if (length(skipped)) do.call(rbind,skipped) |>
  write.csv(file.path(out_dir,'celltypes_not_tested.csv'),row.names=FALSE)
writeLines(c(
  'Unpaired donor-level edgeR quasi-likelihood F test. One library per donor.',
  'Samples from the same donor are not treated as replicates.',
  sprintf('Condition factor levels, sorted: %s', paste(levels, collapse=', ')),
  sprintf('Contrast, coefficient 2: %s', contrast),
  sprintf('A positive logFC means higher in %s than in %s.', levels[2], levels[1]),
  'Cell-type labels are marker-based and unvalidated. An FDR is not confirmed cell-type biology.',
  'FDR_0_05 = 0 means the test ran and no gene passed FDR 0.05.',
  'celltypes_not_tested.csv lists labels that did not meet the donor or count filters.'
), file.path(out_dir,'DESIGN.txt'))
