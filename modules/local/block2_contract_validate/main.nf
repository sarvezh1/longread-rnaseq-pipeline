process BLOCK2_CONTRACT_VALIDATE {
    tag "${meta.id} (Block 2 input contract)"
    label 'process_single'

    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'

    publishDir "${params.outdir}/provenance/samples",
        mode: 'copy',
        pattern: '*.block2_input_validation.json'

    input:
    tuple val(meta), path(bam), path(bai), path(alignment_qc), path(block1_manifest)
    path reference_validation
    path fasta
    path gtf
    val genome_build
    val annotation_release

    output:
    tuple val(meta), path(bam), path(bai), path(alignment_qc), path(block1_manifest), path("${meta.id}.block2_input_validation.json"), emit: validated

    script:
    """
    python ${projectDir}/bin/validate_block2_contract.py \
        --manifest ${block1_manifest} \
        --reference-validation ${reference_validation} \
        --fasta ${fasta} \
        --gtf ${gtf} \
        --sample '${meta.id}' \
        --platform '${meta.platform}' \
        --library-type '${meta.library_type}' \
        --genome-build '${genome_build}' \
        --annotation-release '${annotation_release}' \
        --output ${meta.id}.block2_input_validation.json
    """

    stub:
    """
    printf '{"status":"STUB","sample":"${meta.id}"}\n' > ${meta.id}.block2_input_validation.json
    """
}
