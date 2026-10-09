include { PRIMARY_ISOFORM_DISCOVERY } from '../subworkflows/local/primary_isoform_discovery/main'

workflow BLOCK2_TRANSCRIPTOME {
    take:
    normalized_alignments
    reference_validation
    fasta
    gtf
    genome_build
    annotation_release

    main:
    PRIMARY_ISOFORM_DISCOVERY(
        normalized_alignments,
        reference_validation,
        fasta,
        gtf,
        genome_build,
        annotation_release
    )

    emit:
    transcriptome = PRIMARY_ISOFORM_DISCOVERY.out.transcriptome
    summaries = PRIMARY_ISOFORM_DISCOVERY.out.summaries
    isoquant = PRIMARY_ISOFORM_DISCOVERY.out.isoquant
    sqanti3 = PRIMARY_ISOFORM_DISCOVERY.out.sqanti3
    versions = PRIMARY_ISOFORM_DISCOVERY.out.versions
}
