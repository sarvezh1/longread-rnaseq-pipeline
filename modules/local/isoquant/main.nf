process ISOQUANT {
    tag "${meta.id} (${meta.platform})"
    label 'process_high'

    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'

    publishDir { "${params.outdir}/${meta.benchmark_mode ? 'benchmark/isoquant' : 'isoforms/isoquant/per_sample'}" },
        mode: 'copy',
        pattern: '*.isoquant'
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: '*.isoquant.versions.yml'

    input:
    tuple val(meta), path(bam), path(bai), path(alignment_qc), path(block1_manifest), path(input_validation)
    path fasta
    path gtf

    output:
    tuple val(meta), path("${meta.id}.isoquant"), path("${meta.id}.primary_transcripts.gtf"), path("${meta.id}.transcript_support.tsv.gz"), path("${meta.id}.discovered_transcript_counts.tsv"), path(block1_manifest), path(input_validation), emit: results
    tuple val(meta), path("${meta.id}.isoquant.versions.yml"), emit: versions

    script:
    def data_type = meta.platform == 'ont' ? 'nanopore' : 'pacbio_ccs'
    def full_length = meta.platform == 'pacbio' || meta.preprocessing == 'pychopper'
    def full_length_option = full_length ? '--fl_data' : ''
    def novel_unspliced_option = params.isoquant_report_novel_unspliced == null ? '' : "--report_novel_unspliced ${params.isoquant_report_novel_unspliced}"
    def novel_unspliced_policy = params.isoquant_report_novel_unspliced == null ? 'tool-native' : params.isoquant_report_novel_unspliced.toString()
    """
    set -euo pipefail

    isoquant \
        --reference ${fasta} \
        --genedb ${gtf} \
        --complete_genedb \
        --bam ${bam} \
        --data_type ${data_type} \
        --mode bulk \
        --analysis quantification transcript_discovery exon_quantification \
        --read_group none \
        --large_output read_info read2transcripts \
        --sqanti_output \
        --prefix ${meta.id} \
        --labels ${meta.id} \
        --output ${meta.id}.isoquant \
        --threads ${task.cpus} \
        ${full_length_option} \
        ${novel_unspliced_option}

    model_gtf=\$(find ${meta.id}.isoquant -type f -name '${meta.id}.transcript_models.gtf' -print -quit)
    support=\$(find ${meta.id}.isoquant -type f -name '${meta.id}.transcript_model_reads.tsv.gz' -print -quit)
    counts=\$(find ${meta.id}.isoquant -type f -name '${meta.id}.discovered_transcript_counts.tsv' -print -quit)
    test -n "\${model_gtf}" && test -f "\${model_gtf}"
    test -n "\${support}" && test -s "\${support}"
    test -n "\${counts}" && test -s "\${counts}"
    cp "\${model_gtf}" ${meta.id}.primary_transcripts.gtf
    cp "\${support}" ${meta.id}.transcript_support.tsv.gz
    cp "\${counts}" ${meta.id}.discovered_transcript_counts.tsv

    isoquant_version=\$(isoquant --version 2>&1 | awk 'NF { version=\$NF } END { print version }')
    {
    printf '%s\\n' \
        'ISOQUANT:'
    printf '  isoquant: "%s"\\n' "\${isoquant_version}"
    printf '%s\\n' \
        '  container: "quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142"' \
        '  data_type: "${data_type}"' \
        '  mode: "bulk"' \
        '  analysis: "quantification transcript_discovery exon_quantification"' \
        '  full_length_input: ${full_length}' \
        '  novel_unspliced_policy: "${novel_unspliced_policy}"' \
        '  model_construction_strategy: "tool-native for data type"' \
        '  canonical_junction_annotation: false' \
        '  explicit_minimum_transcript_support: null'
    } > ${meta.id}.isoquant.versions.yml
    """

    stub:
    def data_type = meta.platform == 'ont' ? 'nanopore' : 'pacbio_ccs'
    """
    mkdir -p ${meta.id}.isoquant/${meta.id}
    printf '%s\n' \
        'chrStub\tIsoQuant\ttranscript\t1\t100\t.\t+\t.\tgene_id "stub_gene"; transcript_id "stub_transcript";' \
        'chrStub\tIsoQuant\texon\t1\t100\t.\t+\t.\tgene_id "stub_gene"; transcript_id "stub_transcript"; exon_number "1";' \
        > ${meta.id}.primary_transcripts.gtf
    cp ${meta.id}.primary_transcripts.gtf ${meta.id}.isoquant/${meta.id}/${meta.id}.transcript_models.gtf
    printf 'read_id\ttranscript_id\n' | gzip -c > ${meta.id}.transcript_support.tsv.gz
    printf '#feature_id\tcount\nstub_transcript\t1\n' > ${meta.id}.discovered_transcript_counts.tsv
    cp ${meta.id}.transcript_support.tsv.gz ${meta.id}.isoquant/${meta.id}/${meta.id}.transcript_model_reads.tsv.gz
    cp ${meta.id}.discovered_transcript_counts.tsv ${meta.id}.isoquant/${meta.id}/${meta.id}.discovered_transcript_counts.tsv
    printf '%s\n' \
        'ISOQUANT:' \
        '  isoquant: "4.0.0"' \
        '  container: "quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142"' \
        '  data_type: "${data_type}"' \
        '  mode: "bulk"' \
        > ${meta.id}.isoquant.versions.yml
    """
}
