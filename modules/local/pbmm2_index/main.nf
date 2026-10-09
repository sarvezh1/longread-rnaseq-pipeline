process PBMM2_INDEX {
    tag 'PacBio ISOSEQ index'
    label 'process_high'

    container 'quay.io/biocontainers/pbmm2:26.2.0--h9ee0642_0'

    publishDir "${params.outdir}/alignments/index",
        mode: params.large_file_publish_mode,
        pattern: 'reference.pacbio.isoseq.mmi'
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: 'reference.pacbio.pbmm2_index.versions.yml'

    input:
    path fasta

    output:
    path 'reference.pacbio.isoseq.mmi', emit: index
    path 'reference.pacbio.pbmm2_index.versions.yml', emit: versions

    script:
    """
    pbmm2 index \
        ${fasta} \
        reference.pacbio.isoseq.mmi \
        --preset ISOSEQ \
        --num-threads ${task.cpus}

    cat > reference.pacbio.pbmm2_index.versions.yml <<VERSIONS
    PBMM2_INDEX:
      pbmm2: "\$(pbmm2 --version 2>&1 | sed -n '1s/^pbmm2 //p' | awk '{print \$1}')"
      preset: "ISOSEQ"
      container: "quay.io/biocontainers/pbmm2:26.2.0--h9ee0642_0"
    VERSIONS
    """

    stub:
    """
    touch reference.pacbio.isoseq.mmi
    cat > reference.pacbio.pbmm2_index.versions.yml <<VERSIONS
    PBMM2_INDEX:
      pbmm2: "26.2.0"
      preset: "ISOSEQ"
      container: "quay.io/biocontainers/pbmm2:26.2.0--h9ee0642_0"
    VERSIONS
    """
}

