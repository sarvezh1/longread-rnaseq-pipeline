process ISOSEQ_CLUSTER2 {
    tag "${meta.id} (cluster2)"
    label 'process_high'

    container 'quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0'

    publishDir "${params.outdir}/preprocessing/pacbio",
        mode: params.large_file_publish_mode,
        pattern: "*.clustered.bam*"
    publishDir "${params.outdir}/preprocessing/pacbio",
        mode: 'copy',
        pattern: "*.cluster2.log"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.isoseq_cluster2.versions.yml"

    input:
    tuple val(meta), path(flnc_bam)

    output:
    tuple val(meta), path("${meta.id}.clustered.bam"), emit: clustered_bam
    tuple val(meta), path("${meta.id}.clustered.bam.pbi"), optional: true, emit: pbi
    tuple val(meta), path("${meta.id}.cluster2.log"), emit: log
    tuple val(meta), path("${meta.id}.isoseq_cluster2.versions.yml"), emit: versions

    script:
    def singletons = params.isoseq_include_singletons ? '--singletons' : ''
    """
    isoseq cluster2 \
        ${flnc_bam} \
        ${meta.id}.clustered.bam \
        ${singletons} \
        --log-level INFO \
        > ${meta.id}.cluster2.log 2>&1

    cat > ${meta.id}.isoseq_cluster2.versions.yml <<VERSIONS
    ISOSEQ_CLUSTER2:
      isoseq: "\$(isoseq --version 2>&1 | sed -n '1s/^isoseq //p' | awk '{print \$1}')"
      include_singletons: ${params.isoseq_include_singletons}
      container: "quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0"
    VERSIONS
    """

    stub:
    """
    touch ${meta.id}.clustered.bam
    touch ${meta.id}.clustered.bam.pbi
    touch ${meta.id}.cluster2.log
    cat > ${meta.id}.isoseq_cluster2.versions.yml <<VERSIONS
    ISOSEQ_CLUSTER2:
      isoseq: "26.2.0"
      include_singletons: ${params.isoseq_include_singletons}
      container: "quay.io/biocontainers/isoseq:26.2.0--h9ee0642_0"
    VERSIONS
    """
}
