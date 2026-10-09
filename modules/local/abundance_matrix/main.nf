process ABUNDANCE_MATRIX {
    tag 'production IsoQuant matrices'
    label 'process_single'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/biology", mode: 'copy', pattern: 'abundance'
    input:
    path bundles
    val numerator
    val denominator
    output:
    path 'abundance', emit: dataset
    script:
    """
    python ${projectDir}/bin/build_abundance_matrices.py --bundles ${bundles.join(' ')} --contrast-numerator '${numerator}' --contrast-denominator '${denominator}' --minimum-replicates ${params.minimum_biological_replicates} --design-formula '${params.block4_design_formula}' --output abundance
    """
    stub:
    """
    mkdir abundance
    printf 'feature_id\n' > abundance/gene_counts.tsv; cp abundance/gene_counts.tsv abundance/transcript_counts.tsv
    printf 'feature_id\n' > abundance/gene_tpm.tsv; cp abundance/gene_tpm.tsv abundance/transcript_tpm.tsv
    printf 'sample\tcondition\treplicate\tplatform\tlibrary_type\n' > abundance/sample_metadata.tsv
    printf 'transcript_id\tgene_id\tchrom\tstrand\ttranscript_start\ttranscript_end\ttranscript_length\texon_count\texon_class\texon_coordinates\tintron_chain\tstructural_category\tassociated_gene\tassociated_transcript\tknown_novel\n' > abundance/transcript_metadata.tsv
    printf 'transcript_id\n' > abundance/transcript_presence.tsv; printf 'sample\tcaller_transcript_id\tstandard_transcript_id\tstructure_key\n' > abundance/transcript_id_map.tsv; printf '# empty\n' > abundance/production_transcripts.gtf; : > abundance/production_transcripts.fasta
    printf '{"status":"NOT_RUN_NO_EXPLICIT_CONTRAST","reasons":[]}' > abundance/design_status.json
    """
}
