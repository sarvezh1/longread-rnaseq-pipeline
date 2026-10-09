process TRANSCRIPTOME_CONTRACT {
    tag "${meta.id} (transcriptome contract)"
    label 'process_single'

    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'

    publishDir "${params.outdir}/transcriptome/raw",
        mode: 'copy',
        pattern: '*.raw.{gtf,fasta}'
    publishDir "${params.outdir}/qc/transcript_structure",
        mode: 'copy',
        pattern: '*.{transcript_structure.tsv,transcript_summary.json,transcript_summary.md}'
    publishDir "${params.outdir}/provenance/samples",
        mode: 'copy',
        pattern: '*.transcriptome_contract.json'

    input:
    tuple val(meta), path(isoquant_dir), path(primary_gtf), path(support), path(counts), path(sqanti_dir), path(classification), path(junctions), path(transcript_fasta), path(block1_manifest), path(input_validation), path(isoquant_versions), path(sqanti3_versions)
    path reference_validation
    val genome_build
    val annotation_release

    output:
    tuple val(meta), path("${meta.id}.raw.gtf"), path("${meta.id}.raw.fasta"), path(classification), path(support), path(sqanti_dir), path("${meta.id}.transcript_summary.json"), path("${meta.id}.transcript_structure.tsv"), path("${meta.id}.transcriptome_contract.json"), emit: transcriptome
    tuple val(meta), path("${meta.id}.transcript_summary.md"), emit: human_summary

    script:
    def data_type = meta.platform == 'ont' ? 'nanopore' : 'pacbio_ccs'
    def full_length = meta.platform == 'pacbio' || meta.preprocessing == 'pychopper'
    def novel_unspliced_policy = params.isoquant_report_novel_unspliced == null ? 'tool-native' : params.isoquant_report_novel_unspliced.toString()
    """
    set -euo pipefail
    cp ${primary_gtf} ${meta.id}.raw.gtf
    cp ${transcript_fasta} ${meta.id}.raw.fasta

    python ${projectDir}/bin/summarize_transcript_structure.py \
        --gtf ${primary_gtf} \
        --classification ${classification} \
        --support-counts ${counts} \
        --sample '${meta.id}' \
        --structure-output ${meta.id}.transcript_structure.tsv \
        --json-output ${meta.id}.transcript_summary.json \
        --markdown-output ${meta.id}.transcript_summary.md

    fasta_sha=\$(awk -F'\t' '\$1 == "fasta_sha256" {print \$2}' ${reference_validation})
    annotation_sha=\$(awk -F'\t' '\$1 == "annotation_sha256" {print \$2}' ${reference_validation})
    test -n "\${fasta_sha}" && test -n "\${annotation_sha}"

    python ${projectDir}/bin/build_transcriptome_contract.py \
        --sample '${meta.id}' \
        --platform '${meta.platform}' \
        --library-type '${meta.library_type}' \
        --condition '${meta.condition}' \
        --replicate '${meta.replicate}' \
        --genome-build '${genome_build}' \
        --annotation-release '${annotation_release}' \
        --fasta-sha256 "\${fasta_sha}" \
        --annotation-sha256 "\${annotation_sha}" \
        --isoquant-data-type '${data_type}' \
        --isoquant-full-length '${full_length}' \
        --isoquant-novel-unspliced '${novel_unspliced_policy}' \
        --primary-gtf 'transcriptome/raw/${meta.id}.raw.gtf' \
        --transcript-fasta 'transcriptome/raw/${meta.id}.raw.fasta' \
        --classification 'isoforms/sqanti3/per_sample/${sqanti_dir.name}/${classification.name}' \
        --junctions 'isoforms/sqanti3/per_sample/${sqanti_dir.name}/${junctions.name}' \
        --support-table 'isoforms/isoquant/per_sample/${isoquant_dir.name}/${meta.id}/${meta.id}.transcript_model_reads.tsv.gz' \
        --counts-directory 'isoforms/isoquant/per_sample/${isoquant_dir.name}' \
        --sqanti-directory 'isoforms/sqanti3/per_sample/${sqanti_dir.name}' \
        --structure-table 'qc/transcript_structure/${meta.id}.transcript_structure.tsv' \
        --summary-json 'qc/transcript_structure/${meta.id}.transcript_summary.json' \
        --block1-manifest 'provenance/samples/${block1_manifest.name}' \
        --isoquant-versions 'provenance/software_versions/${isoquant_versions.name}' \
        --sqanti3-versions 'provenance/software_versions/${sqanti3_versions.name}' \
        --output ${meta.id}.transcriptome_contract.json
    """

    stub:
    """
    cp ${primary_gtf} ${meta.id}.raw.gtf
    cp ${transcript_fasta} ${meta.id}.raw.fasta
    printf 'transcript_id\tgene_id\tchrom\tstrand\ttranscript_start\ttranscript_end\ttranscript_length\texon_count\texon_class\texon_coordinates\tintron_chain\tstructural_category\tassociated_gene\tassociated_transcript\tsupporting_read_count\n' > ${meta.id}.transcript_structure.tsv
    printf '{"sample":"${meta.id}","total_transcript_models":1,"single_exon_transcripts":1}\n' > ${meta.id}.transcript_summary.json
    printf '# Transcript structural summary: ${meta.id}\n' > ${meta.id}.transcript_summary.md
    printf '{"contract_version":"2.0","sample":"${meta.id}"}\n' > ${meta.id}.transcriptome_contract.json
    """
}
