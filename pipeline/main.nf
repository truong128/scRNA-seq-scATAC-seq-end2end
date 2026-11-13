nextflow.enable.dsl = 2



process DOWNLOAD_GEO {
  tag "$sid"
  container params.analysis_image
  publishDir "${params.outdir}/inputs", mode: 'copy', overwrite: true
  input:
    tuple val(sid), val(f_url), val(r_url), val(t_url)
    path downloader
  output:
    tuple val(sid), path("${sid}_filtered.h5"), path("${sid}_raw.h5"), path("${sid}_tcr.csv.gz")
  script:
  """
  python ${downloader} --url '${f_url}' --output ${sid}_filtered.h5 --minimum-bytes 10000
  python ${downloader} --url '${r_url}' --output ${sid}_raw.h5 --minimum-bytes 10000
  python ${downloader} --url '${t_url}' --output ${sid}_tcr.csv.gz --minimum-bytes 200
  python - <<'PY'
import gzip, h5py
for name in ['${sid}_filtered.h5', '${sid}_raw.h5']:
    with h5py.File(name, 'r') as f:
        assert 'matrix' in f, 'Not a 10x H5 matrix: '+name
with gzip.open('${sid}_tcr.csv.gz','rt') as f:
    assert 'barcode' in f.readline(), 'Incorrect TCR contig header'
PY
  """
}

process PREP_FASTQ {
  container params.analysis_image
  cpus 2
  memory '4 GB'
  input:
    path prep_code
    path downloader
  output:
    path 'fastqs'
    path 'reference'
  script:
  """
  python ${prep_code} --fastq-url '${params.fastq_url}' --ref-url '${params.ref_url}'
  """
}

process CELLRANGER_COUNT {
  tag 'PBMC1K_CellRanger'

  cpus 8
  memory '48 GB'
  publishDir "${params.outdir}/cellranger", mode: 'copy', pattern: '*.h5'
  input:
    path fq_dir
    path ref_dir
  output:
    tuple val('PBMC1K'), path('PBMC1K_filtered.h5'), path('PBMC1K_raw.h5'), path('PBMC1K_tcr.csv.gz')
  script:
  """
  REF=\$(find \$PWD/reference -name reference.json -print -quit)
  test -n "\$REF"
  REF=\$(dirname "\$REF")
  FQ=\$(find \$PWD/fastqs -name '*R1*fastq.gz' -print -quit)
  test -n "\$FQ"
  FQ=\$(dirname "\$FQ")
  cellranger count --id=PBMC1K_count --transcriptome="\$REF" --fastqs="\$FQ" \
     --sample=pbmc_1k_v3 --localcores=8 --localmem=44 --create-bam=false
  cp PBMC1K_count/outs/filtered_feature_bc_matrix.h5 PBMC1K_filtered.h5
  cp PBMC1K_count/outs/raw_feature_bc_matrix.h5 PBMC1K_raw.h5
  printf 'barcode,chain,productive,raw_clonotype_id\n' | gzip -c > PBMC1K_tcr.csv.gz
  """
}

process ANALYZE {
  tag params.dataset
  container params.analysis_image
  publishDir "${params.outdir}", mode: 'copy', overwrite: true
  cpus 8
  memory '56 GB'
  input:
    path matrices
    path sample_sheet
    path main_code
    path atac_code
    path gene_annotation
    path analysis_yaml
    path marker_csv
    path lr_csv
  output:
    path 'analysis'
  script:
  """
  test -f ${atac_code}
  python ${main_code} --samples ${sample_sheet} --matrix-dir . --config ${analysis_yaml} \
     --markers ${marker_csv} --ligand-receptor ${lr_csv} --outdir analysis \
     --gene-annotation ${gene_annotation} \
     --min-genes ${params.qc_min_genes} --max-genes ${params.qc_max_genes} \
     --max-mt ${params.qc_max_mt} --n-pcs ${params.n_pcs} \
     --resolution ${params.leiden_resolution} \
     --run-scrublet ${params.run_scrublet}
  """
}


process EDGER {
  tag 'replicated_cohort_edgeR'
  publishDir "${params.outdir}/edgeR", mode: 'copy', overwrite: true
  input:
    path analysis_dir
    path rscript
  output:
    path 'edger_results'
  script:
  """
  mkdir -p edger_results
  Rscript ${rscript} ${analysis_dir}/tables edger_results
  """
}

workflow {
  if (params.dataset == 'geo_pilot') {
    records = Channel.fromPath(params.samples).splitCsv(header: true)
      .map { row -> tuple(row.sample_id.toString(), row.filtered_url.toString(),
          row.raw_url.toString(), row.tcr_url.toString()) }
    downloaded = DOWNLOAD_GEO(records, file("${projectDir}/bin/download.py"))
    to_analyze = downloaded.flatMap { sid, filtered, raw, tcr -> [filtered, raw, tcr] }.collect()
    ANALYZE(to_analyze, file(params.samples), file("${projectDir}/bin/analyze.py"),
      file("${projectDir}/bin/atac.py"), file(params.gene_annotation),
      file("${projectDir}/conf/analysis.yaml"),
      file("${projectDir}/assets/cell_markers.csv"), file("${projectDir}/assets/ligand_receptor.csv"))
  } else if (params.dataset == 'pbmc_fastq') {
    if (!params.cellranger_image) {
      error 'pbmc_fastq requires --cellranger_image YOUR_EXISTING_CELLRANGER_IMAGE'
    }
    prepped = PREP_FASTQ(file("${projectDir}/bin/prep_fastq.py"), file("${projectDir}/bin/download.py"))
    counted = CELLRANGER_COUNT(prepped[0], prepped[1])
    ready = counted.flatMap { sid, filtered, raw, tcr -> [filtered, raw, tcr] }.collect()
    ANALYZE(ready, file("${projectDir}/data/fastq_pbmc_1k.csv"),
      file("${projectDir}/bin/analyze.py"), file("${projectDir}/bin/atac.py"), file(params.gene_annotation),
      file("${projectDir}/conf/analysis.yaml"),
      file("${projectDir}/assets/cell_markers.csv"), file("${projectDir}/assets/ligand_receptor.csv"))
  } else if (params.dataset == 'existing') {
    if (!params.local_matrices) {
      error 'existing mode requires --local_matrices /absolute/path/to/h5_files'
    }
    ready = Channel.fromPath("${params.local_matrices}/*_filtered.h5")
      .flatMap { f ->
        def name = f.getName().replace('_filtered.h5','')
        def raw = file("${params.local_matrices}/${name}_raw.h5")
        def tcr = file("${params.local_matrices}/${name}_tcr.csv.gz")
        def peaks = file("${params.local_matrices}/${name}_atac_peaks.h5")
        def single = file("${params.local_matrices}/${name}_atac_singlecell.csv.gz")
        if (!raw.exists()) error "Missing raw 10x H5: ${raw}"
        if (!tcr.exists()) error "Missing TCR placeholder or actual contigs: ${tcr}"
        def out = [f, raw, tcr]
        if (peaks.exists()) out << peaks
        if (single.exists()) out << single
        return out
      }.collect()
    ANALYZE(ready, file(params.samples), file("${projectDir}/bin/analyze.py"),
      file("${projectDir}/bin/atac.py"), file(params.gene_annotation),
      file("${projectDir}/conf/analysis.yaml"), file("${projectDir}/assets/cell_markers.csv"),
      file("${projectDir}/assets/ligand_receptor.csv"))
  } else {
    error "Unknown --dataset ${params.dataset}; geo_pilot, pbmc_fastq, existing"
  }
  if (params.run_edger) {
    EDGER(ANALYZE.out, file("${projectDir}/bin/edger_pseudobulk.R"))
  }
}
