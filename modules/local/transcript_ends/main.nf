process TRANSCRIPT_ENDS {
    tag 'observed transcript boundaries'
    label 'process_single'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'transcript_ends'
    input: path dataset
    output: path 'transcript_ends', emit: results
    script: """python ${projectDir}/bin/transcript_ends.py --dataset ${dataset} --window ${params.transcript_end_cluster_window} --output transcript_ends"""
    stub: """mkdir transcript_ends; printf 'transcript_id\tgene_id\tchrom\tstrand\tobserved_transcript_start\tobserved_transcript_end\tputative_tss\tputative_tes\n' > transcript_ends/transcript_boundaries.tsv; printf 'gene_id\ttranscripts\n' > transcript_ends/gene_boundary_summary.tsv; printf '{"status":"COMPLETE","cluster_window_bp":50}' > transcript_ends/transcript_ends_status.json"""
}
