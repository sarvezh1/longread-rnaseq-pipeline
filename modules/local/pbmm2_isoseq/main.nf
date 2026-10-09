process PBMM2_ISOSEQ {
    tag "${meta.id} (PacBio ISOSEQ)"
    label 'process_high'

    container 'quay.io/biocontainers/pbmm2:26.2.0--h9ee0642_0'

    publishDir "${params.outdir}/alignments/pacbio",
        mode: params.large_file_publish_mode,
        pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/alignments/pacbio",
        mode: 'copy',
        pattern: "*.pbmm2.log"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.pbmm2_isoseq.versions.yml"

    input:
    tuple val(meta), path(clustered_bam)
    path index

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: alignment
    tuple val(meta), path("${meta.id}.pbmm2.log"), emit: log
    tuple val(meta), path("${meta.id}.pbmm2_isoseq.versions.yml"), emit: versions

    script:
    """
    pbmm2 align \
        ${index} \
        ${clustered_bam} \
        ${meta.id}.sorted.bam \
        --preset ISOSEQ \
        --unmapped \
        --sort \
        --bam-index BAI \
        --num-threads ${task.cpus} \
        --sample '${meta.id}' \
        > ${meta.id}.pbmm2.log 2>&1

    cat > ${meta.id}.pbmm2_isoseq.versions.yml <<VERSIONS
    PBMM2_ISOSEQ:
      pbmm2: "\$(pbmm2 --version 2>&1 | sed -n '1s/^pbmm2 //p' | awk '{print \$1}')"
      preset: "ISOSEQ"
      unmapped_records_retained: true
      sample_name: "${meta.id}"
      annotation_guided_alignment: false
      container: "quay.io/biocontainers/pbmm2:26.2.0--h9ee0642_0"
    VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch ${meta.id}.sorted.bam.bai
    touch ${meta.id}.pbmm2.log
    cat > ${meta.id}.pbmm2_isoseq.versions.yml <<VERSIONS
    PBMM2_ISOSEQ:
      pbmm2: "26.2.0"
      preset: "ISOSEQ"
      unmapped_records_retained: true
      sample_name: "${meta.id}"
      annotation_guided_alignment: false
      container: "quay.io/biocontainers/pbmm2:26.2.0--h9ee0642_0"
    VERSIONS
    """
}
