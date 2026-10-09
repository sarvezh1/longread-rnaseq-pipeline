process SATURN_DTU {
    tag 'satuRn differential transcript usage'
    label 'process_medium'
    container 'longread-rnaseq-pipeline/biology:bioc3.23-r4.6'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'transcript_usage'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block4.saturn.versions.yml'
    input: path dataset; val numerator; val denominator; val design_formula
    output: path 'transcript_usage', emit: results; path 'block4.saturn.versions.yml', emit: versions
    script:
    """
    Rscript ${projectDir}/bin/run_saturn.R ${dataset} transcript_usage '${numerator}' '${denominator}' ${params.block4_fdr} '${design_formula}'
    printf 'SATURN:\n  satuRn: "%s"\n  bioconductor: "%s"\n  r: "%s"\n' "\$(Rscript -e 'cat(as.character(packageVersion("satuRn")))')" "\$(Rscript -e 'cat(as.character(BiocManager::version()))')" "\$(Rscript -e 'cat(as.character(getRversion()))')" > block4.saturn.versions.yml
    """
    stub:
    """
    mkdir transcript_usage; printf 'transcript_id\tgene_id\n' > transcript_usage/transcript_usage.tsv; printf 'transcript_id\tgene_id\testimate\tp_value\tadjusted_p_value\tcontrast\n' > transcript_usage/dtu_results.tsv; printf 'transcript_id\tgene_id\ttestable\treason\n' > transcript_usage/dtu_filtering.tsv; printf '{"status":"NOT_RUN_NO_EXPLICIT_CONTRAST"}' > transcript_usage/dtu_status.json; printf '{}' > transcript_usage/versions.json; printf 'SATURN:\n  satuRn: "1.20.0"\n' > block4.saturn.versions.yml
    """
}
