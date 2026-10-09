process BIOLOGY_INTEGRATE {
    tag 'integrated biological evidence'
    label 'process_single'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'integrated'
    input: path dataset; path statistics; path dtu; path splicing; path ends; path orfs; path switches; path enrichment
    output: path 'integrated', emit: results
    script: """python ${projectDir}/bin/integrate_biology.py --dataset ${dataset} --statistics ${statistics} --dtu ${dtu} --splicing ${splicing} --ends ${ends} --orfs ${orfs} --switches ${switches} --enrichment ${enrichment} --output integrated"""
    stub: """mkdir integrated; printf 'gene_id\ttranscript_id\tsample\n' > integrated/integrated_biology.tsv; printf '{"status":"COMPLETE"}' > integrated/integrated_status.json"""
}
