process ISOFORM_SWITCH_ANALYZE {
    tag 'IsoformSwitchAnalyzeR'
    label 'process_high'
    container 'longread-rnaseq-pipeline/biology:bioc3.23-r4.6'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'isoform_switches'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block4.isoform_switch.versions.yml'
    input: path dataset; val numerator; val denominator
    output: path 'isoform_switches', emit: results; path 'block4.isoform_switch.versions.yml', emit: versions
    script:
    """
    Rscript ${projectDir}/bin/run_isoform_switch.R ${dataset} isoform_switches '${numerator}' '${denominator}' ${params.block4_fdr} ${params.isoform_switch_min_dif}
    printf 'ISOFORM_SWITCH_ANALYZER:\n  IsoformSwitchAnalyzeR: "%s"\n  bioconductor: "%s"\n  r: "%s"\n' "\$(Rscript -e 'cat(as.character(packageVersion("IsoformSwitchAnalyzeR")))')" "\$(Rscript -e 'cat(as.character(BiocManager::version()))')" "\$(Rscript -e 'cat(as.character(getRversion()))')" > block4.isoform_switch.versions.yml
    """
    stub:
    """
    mkdir isoform_switches; printf 'gene_id\ttranscript_id\tdIF\tp_value\tadjusted_p_value\tcondition_1\tcondition_2\tstatistical_switch\tfunctional_consequence_evidence\n' > isoform_switches/isoform_switches.tsv; printf '\n' > isoform_switches/predicted_consequences.tsv; printf '{"status":"NOT_RUN_NO_EXPLICIT_CONTRAST"}' > isoform_switches/isoform_switch_status.json; printf '{}' > isoform_switches/versions.json; printf 'ISOFORM_SWITCH_ANALYZER:\n  IsoformSwitchAnalyzeR: "2.12.0"\n' > block4.isoform_switch.versions.yml
    """
}
