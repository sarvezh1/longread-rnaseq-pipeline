process EXPRESSION_QC {
    tag 'expression QC'
    label 'process_low'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'expression_qc'
    input: path dataset
    output: path 'expression_qc', emit: results
    script: """mkdir expression_qc; python ${projectDir}/bin/expression_qc.py --dataset ${dataset} --output expression_qc"""
    stub: """mkdir expression_qc; printf 'sample\tlibrary_count_total\n' > expression_qc/sample_abundance_qc.tsv; printf 'sample\n' > expression_qc/sample_correlation.tsv; printf '{"status":"COMPLETE","pca_status":"NOT_RUN_TOO_FEW_SAMPLES"}' > expression_qc/expression_qc_status.json; : > expression_qc/library_totals.svg; : > expression_qc/detected_transcripts.svg"""
}
