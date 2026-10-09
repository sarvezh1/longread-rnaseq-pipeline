process EDGER_DIFFERENTIAL {
    tag 'edgeR gene/transcript differential expression'
    label 'process_medium'
    container 'longread-rnaseq-pipeline/biology:bioc3.23-r4.6'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'differential_expression'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block4.edger.versions.yml'
    input: path dataset; val numerator; val denominator; val design_formula
    output: path 'differential_expression', emit: results; path 'block4.edger.versions.yml', emit: versions
    script:
    """
    Rscript ${projectDir}/bin/run_edger.R ${dataset} differential_expression '${numerator}' '${denominator}' ${params.block4_fdr} '${design_formula}'
    printf 'EDGER:\n  edgeR: "%s"\n  bioconductor: "%s"\n  r: "%s"\n' "\$(Rscript -e 'cat(as.character(packageVersion("edgeR")))')" "\$(Rscript -e 'cat(as.character(BiocManager::version()))')" "\$(Rscript -e 'cat(as.character(getRversion()))')" > block4.edger.versions.yml
    """
    stub:
    """
    mkdir differential_expression
    printf 'gene_id\tlogFC\tlogCPM\tF\tPValue\tFDR\tcontrast\n' > differential_expression/genes.differential_expression.tsv
    printf 'transcript_id\tlogFC\tlogCPM\tF\tPValue\tFDR\tcontrast\n' > differential_expression/transcripts.differential_expression.tsv
    printf 'feature_id\ttested\treason\n' > differential_expression/genes.filtering.tsv; cp differential_expression/genes.filtering.tsv differential_expression/transcripts.filtering.tsv
    printf '{"status":"NOT_RUN_NO_EXPLICIT_CONTRAST"}' > differential_expression/genes.status.json; cp differential_expression/genes.status.json differential_expression/transcripts.status.json
    printf '{}' > differential_expression/versions.json; printf 'EDGER:\n  edgeR: "4.10.5"\n' > block4.edger.versions.yml
    """
}
