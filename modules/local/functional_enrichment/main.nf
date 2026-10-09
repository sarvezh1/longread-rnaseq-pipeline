process FUNCTIONAL_ENRICHMENT {
    tag 'human GO enrichment'
    label 'process_medium'
    container 'longread-rnaseq-pipeline/biology:bioc3.23-r4.6'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'enrichment'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: 'block4.enrichment.versions.yml'
    input: path dataset; path statistics; path dtu; path switches
    output: path 'enrichment', emit: results; path 'block4.enrichment.versions.yml', emit: versions
    script:
    """
    export BLOCK4_ENABLE_ENRICHMENT='${params.enable_functional_enrichment}' BLOCK4_FDR='${params.block4_fdr}'
    Rscript ${projectDir}/bin/run_enrichment.R ${dataset} ${statistics} ${dtu} ${switches} enrichment
    printf 'ENRICHMENT:\n  clusterProfiler: "%s"\n  org.Hs.eg.db: "%s"\n  bioconductor: "%s"\n  r: "%s"\n' "\$(Rscript -e 'cat(as.character(packageVersion("clusterProfiler")))')" "\$(Rscript -e 'cat(as.character(packageVersion("org.Hs.eg.db")))')" "\$(Rscript -e 'cat(as.character(BiocManager::version()))')" "\$(Rscript -e 'cat(as.character(getRversion()))')" > block4.enrichment.versions.yml
    """
    stub:
    """
    mkdir enrichment; printf 'gene_set\tontology\tID\tDescription\tGeneRatio\tBgRatio\tpvalue\tp.adjust\tqvalue\tgeneID\tCount\n' > enrichment/go_enrichment.tsv; printf 'input_gene_id\tnormalized_ensembl_id\tentrez_id\tmapped\n' > enrichment/identifier_mapping.tsv; printf '{"status":"NOT_RUN_DISABLED"}' > enrichment/enrichment_status.json; printf '{}' > enrichment/versions.json; printf 'ENRICHMENT:\n  clusterProfiler: "4.20.2"\n' > block4.enrichment.versions.yml
    """
}
