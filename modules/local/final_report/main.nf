process FINAL_REPORT {
    tag "${mode} final report"
    label 'process_single'

    container 'ghcr.io/multiqc/multiqc:v1.35@sha256:976be2094a3bc1dab7315ad03440d43aab09bf7756f84f3a92ef7ae1013e4bfd'

    publishDir "${params.outdir}/report", mode: 'copy', pattern: 'final_report/*', saveAs: { name -> name.replaceFirst(/^final_report\//, '') }
    publishDir "${params.outdir}/provenance/final", mode: 'copy', pattern: 'final_provenance/*', saveAs: { name -> name.replaceFirst(/^final_provenance\//, '') }

    input:
    path normalized_contracts, stageAs: 'contracts/block1/*'
    path transcriptome_contracts, stageAs: 'contracts/block2/*'
    path benchmark_contracts, stageAs: 'contracts/block3/*'
    path biology_contracts, stageAs: 'contracts/block4/*'
    path qc_inputs, stageAs: 'qc/source/*'
    path transcript_summaries, stageAs: 'qc/transcript_structure/*'
    path version_files, stageAs: 'versions/*'
    path reference_validation
    path sample_sheet
    path multiqc_dir, stageAs: 'multiqc'
    path status_vocabulary
    path tool_registry
    val parameters_json
    val mode
    val pipeline_version
    val nextflow_version
    val run_id
    val execution_profile
    val started_at

    output:
    path 'final_report/*', emit: report
    path 'final_provenance/*', emit: provenance

    script:
    def parameters_b64 = parameters_json.bytes.encodeBase64().toString()
    """
    printf '%s' '${parameters_b64}' | base64 --decode > resolved_parameters.input.json
    python ${projectDir}/bin/generate_final_report.py \
        --contracts contracts \
        --qc qc \
        --versions versions \
        --reference-validation ${reference_validation} \
        --sample-sheet ${sample_sheet} \
        --parameters resolved_parameters.input.json \
        --status-vocabulary ${status_vocabulary} \
        --tool-registry ${tool_registry} \
        --mode '${mode}' \
        --pipeline-version '${pipeline_version}' \
        --nextflow-version '${nextflow_version}' \
        --run-id '${run_id}' \
        --execution-profile '${execution_profile}' \
        --started-at '${started_at}' \
        --multiqc-report 'multiqc/multiqc_report.html' \
        --output-report final_report \
        --output-provenance final_provenance
    cp -r multiqc final_report/multiqc
    """

    stub:
    """
    mkdir -p final_report/multiqc final_provenance
    printf '<!doctype html><title>Final report stub</title>\n' > final_report/final_report.html
    printf '{"schema_version":"5.0","run":{"overall_status":"COMPLETE"}}\n' > final_report/run_summary.json
    printf 'sample\n' > final_report/samples.tsv
    printf 'analysis\tstatus\tclass\n' > final_report/analysis_status.tsv
    printf 'tool\texpected_version\tobserved_version\tversion_status\n' > final_report/tool_versions.tsv
    cp -r ${multiqc_dir}/* final_report/multiqc/ 2>/dev/null || true
    printf '{"schema_version":"1.0"}\n' > final_provenance/reproducibility_manifest.json
    printf '{"schema_version":"5.0","completion_state":"COMPLETE"}\n' > final_provenance/provenance_bundle.json
    printf '{}\n' > final_provenance/resolved_parameters.json
    printf 'tool\texpected_version\tobserved_version\tversion_status\n' > final_provenance/tool_audit.tsv
    """
}
