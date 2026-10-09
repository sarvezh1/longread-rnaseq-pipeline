process MULTIQC_REPORT {
    tag 'standard QC aggregation'
    label 'process_low'

    container 'ghcr.io/multiqc/multiqc:v1.35@sha256:976be2094a3bc1dab7315ad03440d43aab09bf7756f84f3a92ef7ae1013e4bfd'

    publishDir "${params.outdir}/report", mode: 'copy', pattern: 'multiqc'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block5.multiqc.versions.yml'

    input:
    path qc_inputs, stageAs: 'qc/*'
    path multiqc_config

    output:
    path 'multiqc', emit: report
    path 'block5.multiqc.versions.yml', emit: versions

    script:
    """
    mkdir -p multiqc
    multiqc qc \
        --config ${multiqc_config} \
        --no-ai \
        --force \
        --filename multiqc_report.html \
        --outdir multiqc
    cat > block5.multiqc.versions.yml <<VERSIONS
    MULTIQC:
      MultiQC: "\$(multiqc --version | awk '{print \$NF}')"
      container: "ghcr.io/multiqc/multiqc:v1.35@sha256:976be2094a3bc1dab7315ad03440d43aab09bf7756f84f3a92ef7ae1013e4bfd"
      automated_summaries: "disabled"
      version_check: "disabled"
    VERSIONS
    """

    stub:
    """
    mkdir -p multiqc/multiqc_data
    printf '<!doctype html><title>MultiQC stub</title>\n' > multiqc/multiqc_report.html
    printf '{"report_general_stats_data":[]}\n' > multiqc/multiqc_data/multiqc_data.json
    cat > block5.multiqc.versions.yml <<VERSIONS
    MULTIQC:
      MultiQC: "1.35"
      container: "ghcr.io/multiqc/multiqc:v1.35@sha256:976be2094a3bc1dab7315ad03440d43aab09bf7756f84f3a92ef7ae1013e4bfd"
      automated_summaries: "disabled"
    VERSIONS
    """
}

