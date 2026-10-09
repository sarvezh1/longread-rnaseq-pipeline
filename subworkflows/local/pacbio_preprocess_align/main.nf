include { NANOPLOT_RAW as NANOPLOT_RAW_PACBIO } from '../../../modules/local/nanoplot_raw/main'
include { LIMA_ISOSEQ }                         from '../../../modules/local/lima_isoseq/main'
include { ISOSEQ_REFINE }                       from '../../../modules/local/isoseq_refine/main'
include { ISOSEQ_CLUSTER2 }                     from '../../../modules/local/isoseq_cluster2/main'
include { PBMM2_INDEX }                         from '../../../modules/local/pbmm2_index/main'
include { PBMM2_ISOSEQ }                        from '../../../modules/local/pbmm2_isoseq/main'

workflow PACBIO_PREPROCESS_ALIGN {
    take:
    samples
    fasta

    main:
    NANOPLOT_RAW_PACBIO(samples.map { meta, reads, _primers -> tuple(meta, reads) })

    input_routes = samples.branch { sample_record ->
        hifi: sample_record[0].library_type == 'isoseq_hifi'
        flnc: sample_record[0].library_type == 'isoseq_flnc'
    }

    LIMA_ISOSEQ(input_routes.hifi)
    ISOSEQ_REFINE(LIMA_ISOSEQ.out.full_length_bam)

    hifi_flnc = ISOSEQ_REFINE.out.flnc_bam
    supplied_flnc = input_routes.flnc.map { meta, reads, _primers -> tuple(meta, reads) }
    flnc_reads = hifi_flnc.mix(supplied_flnc)

    ISOSEQ_CLUSTER2(flnc_reads)

    index_trigger = samples
        .map { _meta, _reads, _primers -> true }
        .first()
        .combine(fasta)
        .map { _trigger, reference -> reference }
    PBMM2_INDEX(index_trigger)
    pacbio_index = PBMM2_INDEX.out.index.first()

    PBMM2_ISOSEQ(ISOSEQ_CLUSTER2.out.clustered_bam, pacbio_index)

    versions = NANOPLOT_RAW_PACBIO.out.versions
        .mix(LIMA_ISOSEQ.out.versions)
        .mix(ISOSEQ_REFINE.out.versions)
        .mix(ISOSEQ_CLUSTER2.out.versions)
        .mix(PBMM2_ISOSEQ.out.versions)

    emit:
    alignment = PBMM2_ISOSEQ.out.alignment
    raw_qc = NANOPLOT_RAW_PACBIO.out.qc
    preprocessing_logs = LIMA_ISOSEQ.out.log
        .mix(ISOSEQ_REFINE.out.log)
        .mix(ISOSEQ_CLUSTER2.out.log)
    versions = versions
}
