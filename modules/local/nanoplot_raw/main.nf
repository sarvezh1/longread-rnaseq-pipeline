process NANOPLOT_RAW {
    tag "${meta.id} (${meta.platform})"
    label 'process_low'

    container 'quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_1'

    publishDir "${params.outdir}/qc/raw",
        mode: 'copy',
        pattern: "*.nanoplot"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.nanoplot.versions.yml"

    input:
    tuple val(meta), path(reads)

    output:
    tuple val(meta), path("${meta.id}.nanoplot"), emit: qc
    tuple val(meta), path("${meta.id}.nanoplot.versions.yml"), emit: versions

    script:
    def input_mode = meta.platform == 'ont' ? '--fastq' : '--ubam'
    """
    NanoPlot \
        ${input_mode} ${reads} \
        --threads ${task.cpus} \
        --tsv_stats \
        --info_in_report \
        --outdir ${meta.id}.nanoplot

    cat > ${meta.id}.nanoplot.versions.yml <<VERSIONS
    NANOPLOT_RAW:
      NanoPlot: "\$(NanoPlot --version 2>&1 | cut -d' ' -f2)"
      container: "quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_1"
    VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.nanoplot
    touch ${meta.id}.nanoplot/NanoPlot-report.html
    touch ${meta.id}.nanoplot/NanoStats.txt
    cat > ${meta.id}.nanoplot.versions.yml <<VERSIONS
    NANOPLOT_RAW:
      NanoPlot: "1.48.0"
      container: "quay.io/biocontainers/nanoplot:1.48.0--pyhdfd78af_1"
    VERSIONS
    """
}
