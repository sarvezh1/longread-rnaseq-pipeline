process SAMPLE_ABUNDANCE {
    tag "${meta.id} (production abundance)"
    label 'process_single'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/biology/abundance/per_sample", mode: 'copy', pattern: '*.abundance'
    input:
    tuple val(meta), path(gtf), path(fasta), path(classification), path(support), path(sqanti_dir), path(summary), path(structure), path(contract), path(isoquant_dir)
    output:
    tuple val(meta), path("${meta.id}.abundance"), emit: bundle
    script:
    """
    python ${projectDir}/bin/prepare_sample_abundance.py --sample '${meta.id}' --platform '${meta.platform}' --library-type '${meta.library_type}' --condition '${meta.condition}' --replicate '${meta.replicate}' --isoquant-dir ${isoquant_dir} --gtf ${gtf} --fasta ${fasta} --structure ${structure} --classification ${classification} --contract ${contract} --output ${meta.id}.abundance
    """
    stub:
    """
    mkdir -p ${meta.id}.abundance
    printf '{"sample":"${meta.id}","condition":"${meta.condition}","replicate":"${meta.replicate}","platform":"${meta.platform}","library_type":"${meta.library_type}"}\n' > ${meta.id}.abundance/sample_metadata.json
    printf 'feature_id\tcount\n' > ${meta.id}.abundance/gene_counts.tsv
    printf 'feature_id\tcount\n' > ${meta.id}.abundance/transcript_counts.tsv
    printf 'feature_id\tTPM\n' > ${meta.id}.abundance/gene_tpm.tsv
    printf 'feature_id\tTPM\n' > ${meta.id}.abundance/transcript_tpm.tsv
    cp ${gtf} ${meta.id}.abundance/transcripts.gtf; cp ${fasta} ${meta.id}.abundance/transcripts.fasta; cp ${structure} ${meta.id}.abundance/transcript_structure.tsv; cp ${classification} ${meta.id}.abundance/sqanti3_classification.tsv; cp ${contract} ${meta.id}.abundance/block2_contract.json
    """
}
