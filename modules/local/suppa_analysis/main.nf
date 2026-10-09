process SUPPA_ANALYSIS {
    tag 'SUPPA events and PSI'
    label 'process_medium'
    container 'longread-rnaseq-pipeline/suppa:2.4'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'alternative_splicing'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block4.suppa.versions.yml'
    input: path dataset; val numerator; val denominator
    output: path 'alternative_splicing', emit: results; path 'block4.suppa.versions.yml', emit: versions
    script:
    """
    python ${projectDir}/bin/run_suppa.py --dataset ${dataset} --output suppa_raw --numerator '${numerator}' --denominator '${denominator}'
    mkdir -p alternative_splicing/events alternative_splicing/psi alternative_splicing/differential
    cp suppa_raw/events.tsv suppa_raw/events_status.json alternative_splicing/events/
    cp suppa_raw/psi.tsv suppa_raw/psi_status.json alternative_splicing/psi/
    cp suppa_raw/differential_splicing.tsv suppa_raw/differential_status.json alternative_splicing/differential/
    mv suppa_raw alternative_splicing/native
    printf '{"release":"2.4","cli_reported":"%s","container":"longread-rnaseq-pipeline/suppa:2.4"}\n' "\$(suppa.py --version 2>&1 | awk '{print \$NF}')" > alternative_splicing/runtime_versions.json
    printf 'SUPPA:\n  release: "2.4"\n  cli_reported: "%s"\n  container: "longread-rnaseq-pipeline/suppa:2.4"\n' "\$(suppa.py --version 2>&1 | awk '{print \$NF}')" > block4.suppa.versions.yml
    """
    stub:
    """
    mkdir -p alternative_splicing/events alternative_splicing/psi alternative_splicing/differential alternative_splicing/native
    printf 'event_id\tgene_id\tevent_type\talternative_transcripts\ttotal_transcripts\n' > alternative_splicing/events/events.tsv; printf '{"status":"COMPLETE_EMPTY_RESULT"}' > alternative_splicing/events/events_status.json
    printf 'event_id\tgene_id\tevent_type\talternative_transcripts\ttotal_transcripts\n' > alternative_splicing/psi/psi.tsv; printf '{"status":"NOT_RUN_NO_EVENTS"}' > alternative_splicing/psi/psi_status.json
    printf 'event_id\tgene_id\tevent_type\talternative_transcripts\ttotal_transcripts\tdelta_psi\tp_value\tadjusted_p_value\tcontrast\n' > alternative_splicing/differential/differential_splicing.tsv; cp alternative_splicing/psi/psi_status.json alternative_splicing/differential/differential_status.json; printf '{"release":"2.4","cli_reported":"2.3"}' > alternative_splicing/runtime_versions.json; printf 'SUPPA:\n  release: "2.4"\n' > block4.suppa.versions.yml
    """
}
