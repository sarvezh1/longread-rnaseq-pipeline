process TRANSCRIPT_STANDARDIZE {
    tag "${meta.benchmark_sample_id} (${meta.benchmark_caller})"
    label 'process_single'

    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'

    publishDir "${params.outdir}/benchmark/concordance/transcripts", mode: 'copy', pattern: '*.standardized_transcripts.tsv'
    publishDir "${params.outdir}/benchmark/concordance/splice_junctions", mode: 'copy', pattern: '*.junctions.tsv'
    publishDir "${params.outdir}/benchmark/summary", mode: 'copy', pattern: '*.caller_summary.json'

    input:
    tuple val(meta), path(caller_dir), path(primary_gtf), path(support), path(counts), path(sqanti_dir), path(classification), path(sqanti_junctions), path(transcript_fasta), path(block1_manifest), path(input_validation)

    output:
    tuple val(meta), path("${meta.id}.standardized_transcripts.tsv"), path("${meta.id}.junctions.tsv"), path("${meta.id}.caller_summary.json"), path(block1_manifest), emit: standardized

    script:
    """
    python ${projectDir}/bin/standardize_benchmark_transcripts.py \
        --sample '${meta.benchmark_sample_id}' \
        --platform '${meta.platform}' \
        --condition '${meta.condition}' \
        --replicate '${meta.replicate}' \
        --caller '${meta.benchmark_caller}' \
        --gtf ${primary_gtf} \
        --classification ${classification} \
        --sqanti-junctions ${sqanti_junctions} \
        --support ${support} \
        --support-kind '${meta.benchmark_caller == 'flair' ? 'count' : 'auto'}' \
        --transcripts ${meta.id}.standardized_transcripts.tsv \
        --junctions ${meta.id}.junctions.tsv \
        --summary ${meta.id}.caller_summary.json
    """

    stub:
    """
    printf 'sample\tplatform\tcondition\treplicate\tcaller\tcaller_transcript_id\tchrom\tstrand\ttranscript_start\ttranscript_end\ttss\ttes\ttranscript_length\texon_count\texon_class\texon_coordinates\tintron_chain\tgene_assignment\tcaller_gene_id\tsqanti_category\tsqanti_associated_transcript\tsupport_value\tsupport_semantics\n' > ${meta.id}.standardized_transcripts.tsv
    printf 'sample\tplatform\tcaller\tcaller_transcript_id\tchrom\tstrand\tintron_start\tintron_end\tcanonical_status\treference_status\tsupport_value\tsupport_semantics\n' > ${meta.id}.junctions.tsv
    printf '{"sample":"${meta.benchmark_sample_id}","platform":"${meta.platform}","caller":"${meta.benchmark_caller}","transcripts":0}\n' > ${meta.id}.caller_summary.json
    """
}
