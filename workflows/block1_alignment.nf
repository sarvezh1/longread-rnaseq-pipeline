include { ONT_PREPROCESS_ALIGN }     from '../subworkflows/local/ont_preprocess_align/main'
include { PACBIO_PREPROCESS_ALIGN }  from '../subworkflows/local/pacbio_preprocess_align/main'
include { ALIGNMENT_QC }             from '../subworkflows/local/alignment_qc/main'
include { REFERENCE_VALIDATE }       from '../modules/local/reference_validate/main'

workflow BLOCK1_ALIGNMENT {
    take:
    samples
    fasta
    fasta_source
    gtf_source
    genome_build
    annotation_release

    main:
    REFERENCE_VALIDATE(fasta, gtf_source)

    platform_routes = samples.branch { sample_record ->
        ont: sample_record[0].platform == 'ont'
        pacbio: sample_record[0].platform == 'pacbio'
    }

    ONT_PREPROCESS_ALIGN(platform_routes.ont, REFERENCE_VALIDATE.out.validated_fasta)
    PACBIO_PREPROCESS_ALIGN(platform_routes.pacbio, REFERENCE_VALIDATE.out.validated_fasta)

    aligned_reads = ONT_PREPROCESS_ALIGN.out.alignment
        .mix(PACBIO_PREPROCESS_ALIGN.out.alignment)
    runtime_versions = ONT_PREPROCESS_ALIGN.out.versions
        .mix(PACBIO_PREPROCESS_ALIGN.out.versions)

    ALIGNMENT_QC(
        aligned_reads,
        REFERENCE_VALIDATE.out.validated_fasta,
        REFERENCE_VALIDATE.out.fai,
        runtime_versions,
        REFERENCE_VALIDATE.out.report,
        fasta_source,
        gtf_source,
        genome_build,
        annotation_release
    )

    emit:
    normalized = ALIGNMENT_QC.out.normalized
    manifests = ALIGNMENT_QC.out.manifests
    alignment_qc = ALIGNMENT_QC.out.qc
    ont_raw_qc = ONT_PREPROCESS_ALIGN.out.raw_qc
    pacbio_raw_qc = PACBIO_PREPROCESS_ALIGN.out.raw_qc
    reference_validation = REFERENCE_VALIDATE.out.report
    reference_versions = REFERENCE_VALIDATE.out.versions
    runtime_versions = ALIGNMENT_QC.out.runtime_versions
}
