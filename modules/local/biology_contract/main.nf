process BIOLOGY_CONTRACT {
    tag 'Block 4 contract'
    label 'process_single'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'biology_contract.json'
    input: path dataset; path qc; path statistics; path dtu; path splicing; path ends; path orfs; path switches; path enrichment; path integrated
    output: path 'biology_contract.json', emit: contract
    script: """python ${projectDir}/bin/build_biology_contract.py --dataset ${dataset} --qc ${qc} --statistics ${statistics} --dtu ${dtu} --splicing ${splicing} --ends ${ends} --orfs ${orfs} --switches ${switches} --enrichment ${enrichment} --integrated ${integrated} --output biology_contract.json"""
    stub: """printf '{"contract_version":"4.0","production_transcriptome":"IsoQuant via Block 2 standardized transcriptome contract"}' > biology_contract.json"""
}
