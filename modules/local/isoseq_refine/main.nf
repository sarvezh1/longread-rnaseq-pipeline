process ISOSEQ_REFINE {
    tag "${meta.id} (FLNC)"
    label 'process_high'

    container 'quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0'

    publishDir "${params.outdir}/preprocessing/pacbio",
        mode: params.large_file_publish_mode,
        pattern: "*.flnc.bam*"
    publishDir "${params.outdir}/preprocessing/pacbio",
        mode: 'copy',
        pattern: "*.refine.log"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.isoseq_refine.versions.yml"

    input:
    tuple val(meta), path(full_length_bam), path(primers)

    output:
    tuple val(meta), path("${meta.id}.flnc.bam"), emit: flnc_bam
    tuple val(meta), path("${meta.id}.flnc.bam.pbi"), optional: true, emit: pbi
    tuple val(meta), path("${meta.id}.refine.log"), emit: log
    tuple val(meta), path("${meta.id}.isoseq_refine.versions.yml"), emit: versions

    script:
    def require_polya = meta.require_polya ? '--require-polya' : ''
    """
    isoseq refine \
        ${full_length_bam} \
        ${primers} \
        ${meta.id}.flnc.bam \
        ${require_polya} \
        --log-level INFO \
        > ${meta.id}.refine.log 2>&1

    cat > ${meta.id}.isoseq_refine.versions.yml <<VERSIONS
    ISOSEQ_REFINE:
      isoseq: "\$(isoseq --version 2>&1 | sed -n '1s/^isoseq //p' | awk '{print \$1}')"
      require_polya: ${meta.require_polya}
      container: "quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0"
    VERSIONS
    """

    stub:
    """
    touch ${meta.id}.flnc.bam
    touch ${meta.id}.flnc.bam.pbi
    touch ${meta.id}.refine.log
    cat > ${meta.id}.isoseq_refine.versions.yml <<VERSIONS
    ISOSEQ_REFINE:
      isoseq: "26.2.0"
      require_polya: ${meta.require_polya}
      container: "quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0"
    VERSIONS
    """
}
