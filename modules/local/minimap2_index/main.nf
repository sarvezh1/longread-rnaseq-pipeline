process MINIMAP2_INDEX {
    tag 'ONT splice index'
    label 'process_high'

    container 'community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce'

    publishDir "${params.outdir}/alignments/index",
        mode: params.large_file_publish_mode,
        pattern: 'reference.ont.splice.mmi'
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: 'reference.ont.minimap2_index.versions.yml'

    input:
    path fasta

    output:
    path 'reference.ont.splice.mmi', emit: index
    path 'reference.ont.minimap2_index.versions.yml', emit: versions

    script:
    """
    minimap2 \
        -x splice \
        -t ${task.cpus} \
        -d reference.ont.splice.mmi \
        ${fasta}

    cat > reference.ont.minimap2_index.versions.yml <<VERSIONS
    MINIMAP2_INDEX:
      minimap2: "\$(minimap2 --version)"
      samtools: "\$(samtools version | sed -n '1s/^samtools //p')"
      container: "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    VERSIONS
    """

    stub:
    """
    touch reference.ont.splice.mmi
    cat > reference.ont.minimap2_index.versions.yml <<VERSIONS
    MINIMAP2_INDEX:
      minimap2: "2.30"
      samtools: "1.23.1"
      container: "community.wave.seqera.io/library/minimap2_samtools:b09096fc890429ce"
    VERSIONS
    """
}

