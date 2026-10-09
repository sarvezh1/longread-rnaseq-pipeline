process TRANSCRIPT_CONCORDANCE {
    tag "${meta.id} (${meta.platform})"
    label 'process_single'

    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'

    publishDir "${params.outdir}/benchmark/concordance/membership", mode: 'copy', pattern: '*.membership.tsv'
    publishDir "${params.outdir}/benchmark/concordance/intron_chains", mode: 'copy', pattern: '*.pairwise_metrics.tsv'
    publishDir "${params.outdir}/benchmark/concordance/transcript_ends", mode: 'copy', pattern: '*.transcript_ends.tsv'
    publishDir "${params.outdir}/benchmark/concordance/splice_junctions", mode: 'copy', pattern: '*.partial_junction_overlap.tsv'
    publishDir "${params.outdir}/benchmark/concordance/transcripts", mode: 'copy', pattern: '*.gene_overlap.tsv'
    publishDir "${params.outdir}/benchmark/summary", mode: 'copy', pattern: '*.{benchmark_summary.json,benchmark_contract.json}'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: '*.benchmark.versions.yml'

    input:
    tuple val(meta), path(isoquant_structures), path(flair_structures), path(bambu_structures), path(block1_manifest), path(isoquant_versions), path(flair_versions), path(bambu_versions), path(sqanti_isoquant_versions), path(sqanti_flair_versions), path(sqanti_bambu_versions)

    output:
    tuple val(meta), path("${meta.id}.benchmark_summary.json"), path("${meta.id}.benchmark_contract.json"), path("${meta.id}.membership.tsv"), path("${meta.id}.pairwise_metrics.tsv"), path("${meta.id}.transcript_ends.tsv"), path("${meta.id}.partial_junction_overlap.tsv"), path("${meta.id}.gene_overlap.tsv"), emit: benchmark
    tuple val(meta), path("${meta.id}.benchmark.versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    python ${projectDir}/bin/compare_transcript_structures.py \
        --sample '${meta.id}' --platform '${meta.platform}' --condition '${meta.condition}' --replicate '${meta.replicate}' \
        --isoquant ${isoquant_structures} --flair ${flair_structures} --bambu ${bambu_structures} \
        --end-bins '${params.benchmark_end_bins}' --single-exon-overlap ${params.benchmark_single_exon_reciprocal_overlap} \
        --membership ${meta.id}.membership.tsv \
        --pairwise ${meta.id}.pairwise_metrics.tsv \
        --ends ${meta.id}.transcript_ends.tsv \
        --partial-junctions ${meta.id}.partial_junction_overlap.tsv \
        --gene-overlap ${meta.id}.gene_overlap.tsv \
        --summary ${meta.id}.benchmark_summary.json

    {
      printf 'BENCHMARK_RUNTIME:\n'
      printf '  isoquant:\n    caller:\n'
      sed '1d; s/^/      /' ${isoquant_versions}
      printf '    sqanti3:\n'
      sed '1d; s/^/      /' ${sqanti_isoquant_versions}
      printf '  flair:\n    caller:\n'
      sed '1d; s/^/      /' ${flair_versions}
      printf '    sqanti3:\n'
      sed '1d; s/^/      /' ${sqanti_flair_versions}
      printf '  bambu:\n    caller:\n'
      sed '1d; s/^/      /' ${bambu_versions}
      printf '    sqanti3:\n'
      sed '1d; s/^/      /' ${sqanti_bambu_versions}
    } > ${meta.id}.benchmark.versions.yml

    python ${projectDir}/bin/build_benchmark_contract.py \
        --sample '${meta.id}' --platform '${meta.platform}' --condition '${meta.condition}' --replicate '${meta.replicate}' \
        --summary ${meta.id}.benchmark_summary.json --manifest ${block1_manifest} \
        --version ${isoquant_versions} --version ${flair_versions} --version ${bambu_versions} \
        --version ${sqanti_isoquant_versions} --version ${sqanti_flair_versions} --version ${sqanti_bambu_versions} \
        --output ${meta.id}.benchmark_contract.json
    """

    stub:
    """
    printf '{"sample":"${meta.id}","platform":"${meta.platform}","contract_version":"3.0"}\n' > ${meta.id}.benchmark_summary.json
    printf '{"sample":"${meta.id}","contract_version":"3.0","consensus_transcriptome_created":false}\n' > ${meta.id}.benchmark_contract.json
    printf 'sample\tplatform\tstructure_group\tmatching_basis\tchrom\tstrand\tisoquant_transcript_ids\tflair_transcript_ids\tbambu_transcript_ids\tmembership\tsqanti_category_groups\n' > ${meta.id}.membership.tsv
    printf 'sample\tplatform\tcaller_pair\tlevel\tintersection\tunion\tjaccard\tstatus\n' > ${meta.id}.pairwise_metrics.tsv
    printf 'sample\tplatform\tcaller_pair\n' > ${meta.id}.transcript_ends.tsv
    printf 'sample\tplatform\tcaller_pair\n' > ${meta.id}.partial_junction_overlap.tsv
    printf 'sample\tplatform\tcaller_pair\n' > ${meta.id}.gene_overlap.tsv
    printf 'BENCHMARK_RUNTIME:\n  status: "STUB"\n' > ${meta.id}.benchmark.versions.yml
    """
}
