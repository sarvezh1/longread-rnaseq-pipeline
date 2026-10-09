process PYCHOPPER {
    tag "${meta.id} (${meta.library_type})"
    label 'process_medium'

    container 'quay.io/biocontainers/pychopper:2.7.10--pyhdfd78af_0'

    publishDir "${params.outdir}/preprocessing/ont",
        mode: params.large_file_publish_mode,
        pattern: "*.fastq"
    publishDir "${params.outdir}/preprocessing/ont",
        mode: 'copy',
        pattern: "*.{pdf,tsv}"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.pychopper.versions.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.full_length.fastq"), path("${meta.id}.rescued.fastq"), emit: processed_reads
    tuple val(meta), path("${meta.id}.unclassified.fastq"), path("${meta.id}.quality_fail.fastq"), path("${meta.id}.length_fail.fastq"), emit: rejected_reads
    tuple val(meta), path("${meta.id}.pychopper_report.pdf"), path("${meta.id}.pychopper_stats.tsv"), path("${meta.id}.pychopper_reads.tsv"), emit: reports
    tuple val(meta), path("${meta.id}.pychopper.versions.yml"), emit: versions

    script:
    """
    pychopper \
        -k ${meta.ont_kit} \
        -Q 0 \
        -z 0 \
        -t ${task.cpus} \
        -r ${meta.id}.pychopper_report.pdf \
        -S ${meta.id}.pychopper_stats.tsv \
        -D ${meta.id}.pychopper_reads.tsv \
        -u ${meta.id}.unclassified.fastq \
        -w ${meta.id}.rescued.fastq \
        -K ${meta.id}.quality_fail.fastq \
        -l ${meta.id}.length_fail.fastq \
        ${reads} \
        ${meta.id}.full_length.fastq

    cat > ${meta.id}.pychopper.versions.yml <<VERSIONS
    PYCHOPPER:
      pychopper: "\$(python -c 'from importlib.metadata import version; print(version("pychopper"))')"
      container: "quay.io/biocontainers/pychopper:2.7.10--pyhdfd78af_0"
      minimum_read_quality: 0
      minimum_fragment_length: 0
    VERSIONS
    """

    stub:
    """
    touch ${meta.id}.full_length.fastq
    touch ${meta.id}.rescued.fastq
    touch ${meta.id}.unclassified.fastq
    touch ${meta.id}.quality_fail.fastq
    touch ${meta.id}.length_fail.fastq
    touch ${meta.id}.pychopper_report.pdf
    touch ${meta.id}.pychopper_stats.tsv
    touch ${meta.id}.pychopper_reads.tsv
    cat > ${meta.id}.pychopper.versions.yml <<VERSIONS
    PYCHOPPER:
      pychopper: "2.7.10"
      container: "quay.io/biocontainers/pychopper:2.7.10--pyhdfd78af_0"
      minimum_read_quality: 0
      minimum_fragment_length: 0
    VERSIONS
    """
}
