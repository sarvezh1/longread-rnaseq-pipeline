include { MULTIQC_REPORT } from '../../../modules/local/multiqc_report/main'
include { FINAL_REPORT }   from '../../../modules/local/final_report/main'

workflow FINAL_REPORTING {
    take:
    normalized_contracts
    transcriptome_contracts
    benchmark_contracts
    biology_contracts
    qc_inputs
    transcript_summaries
    version_files
    reference_validation
    sample_sheet
    status_vocabulary
    tool_registry
    multiqc_config
    parameters_json
    mode
    pipeline_version
    nextflow_version
    run_id
    execution_profile
    started_at

    main:
    MULTIQC_REPORT(qc_inputs.collect(), multiqc_config)
    report_versions = version_files.mix(MULTIQC_REPORT.out.versions)
    FINAL_REPORT(
        normalized_contracts.collect(),
        transcriptome_contracts.collect(),
        benchmark_contracts.collect(),
        biology_contracts.collect(),
        qc_inputs.collect(),
        transcript_summaries.collect(),
        report_versions.collect(),
        reference_validation,
        sample_sheet,
        MULTIQC_REPORT.out.report,
        status_vocabulary,
        tool_registry,
        parameters_json,
        mode,
        pipeline_version,
        nextflow_version,
        run_id,
        execution_profile,
        started_at
    )

    emit:
    report = FINAL_REPORT.out.report
    provenance = FINAL_REPORT.out.provenance
    multiqc = MULTIQC_REPORT.out.report
    versions = MULTIQC_REPORT.out.versions
}
