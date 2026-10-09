process MINIMAP2_SPLICE {
    tag "${meta.id} (ONT splice)"
    label 'process_high'

    container 'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce'

    publishDir "${params.outdir}/alignments/ont",
        mode: params.large_file_publish_mode,
        pattern: "*.{bam,bai}"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.minimap2_splice.versions.yml"

    input:
    tuple val(meta), path(reads)
    path index

    output:
    tuple val(meta), path("${meta.id}.sorted.bam"), path("${meta.id}.sorted.bam.bai"), emit: alignment
    tuple val(meta), path("${meta.id}.minimap2_splice.versions.yml"), emit: versions

    script:
    def sort_threads = Math.max(1, (task.cpus as int) - 1)
    """
    set -o pipefail

    minimap2 \
        -ax splice \
        -t ${task.cpus} \
        -R '@RG\\tID:${meta.id}\\tSM:${meta.id}\\tPL:ONT\\tLB:${meta.library_type}' \
        ${index} \
        ${reads} \
    | samtools sort \
        -@ ${sort_threads} \
        -o ${meta.id}.sorted.bam \
        -

    samtools index -@ ${task.cpus} ${meta.id}.sorted.bam

    cat > ${meta.id}.minimap2_splice.versions.yml <<VERSIONS
    MINIMAP2_SPLICE:
      minimap2: "\$(minimap2 --version)"
      samtools: "\$(samtools version | sed -n '1s/^samtools //p')"
      preset: "splice"
      annotation_guided_alignment: false
      container: "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    VERSIONS
    """

    stub:
    """
    touch ${meta.id}.sorted.bam
    touch ${meta.id}.sorted.bam.bai
    cat > ${meta.id}.minimap2_splice.versions.yml <<VERSIONS
    MINIMAP2_SPLICE:
      minimap2: "2.30"
      samtools: "1.23.1"
      preset: "splice"
      annotation_guided_alignment: false
      container: "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    VERSIONS
    """
}
