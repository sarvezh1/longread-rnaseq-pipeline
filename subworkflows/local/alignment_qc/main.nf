include { SAMTOOLS_ALIGNMENT_QC } from '../../../modules/local/samtools_alignment_qc/main'
include { NORMALIZED_MANIFEST }   from '../../../modules/local/normalized_manifest/main'

workflow ALIGNMENT_QC {
    take:
    alignments
    fasta
    fasta_index
    runtime_versions
    reference_validation
    fasta_source
    gtf_source
    genome_build
    annotation_release

    main:
    SAMTOOLS_ALIGNMENT_QC(alignments, fasta, fasta_index)

    keyed_alignments = alignments.map { meta, bam, bai ->
        tuple(meta.id, meta, bam, bai)
    }
    keyed_qc = SAMTOOLS_ALIGNMENT_QC.out.qc.map { meta, qc_dir ->
        tuple(meta.id, qc_dir)
    }
    aligned_with_qc = keyed_alignments
        .join(keyed_qc)
        .map { _sample_id, meta, bam, bai, qc_dir -> tuple(meta, bam, bai, qc_dir) }

    all_runtime_versions = runtime_versions.mix(SAMTOOLS_ALIGNMENT_QC.out.versions)
    keyed_versions = all_runtime_versions
        .map { meta, versions_file -> tuple(meta.id, versions_file) }
        .groupTuple()
    normalized_inputs = aligned_with_qc
        .map { meta, bam, bai, qc_dir -> tuple(meta.id, meta, bam, bai, qc_dir) }
        .join(keyed_versions)
        .map { _sample_id, meta, bam, bai, qc_dir, versions_files -> tuple(meta, bam, bai, qc_dir, versions_files) }

    NORMALIZED_MANIFEST(
        normalized_inputs,
        reference_validation,
        fasta_source,
        gtf_source,
        genome_build,
        annotation_release
    )

    emit:
    normalized = NORMALIZED_MANIFEST.out.normalized
    manifests = NORMALIZED_MANIFEST.out.manifest_tsv
    qc = SAMTOOLS_ALIGNMENT_QC.out.qc
    versions = SAMTOOLS_ALIGNMENT_QC.out.versions
    runtime_versions = NORMALIZED_MANIFEST.out.runtime_versions
}
