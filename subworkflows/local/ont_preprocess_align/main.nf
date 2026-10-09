include { NANOPLOT_RAW as NANOPLOT_RAW_ONT } from '../../../modules/local/nanoplot_raw/main'
include { PYCHOPPER }                          from '../../../modules/local/pychopper/main'
include { MINIMAP2_INDEX }                     from '../../../modules/local/minimap2_index/main'
include { MINIMAP2_SPLICE }                    from '../../../modules/local/minimap2_splice/main'

workflow ONT_PREPROCESS_ALIGN {
    take:
    samples
    fasta

    main:
    NANOPLOT_RAW_ONT(samples.map { meta, reads, _primers -> tuple(meta, reads) })

    preprocessing_routes = samples.branch { sample_record ->
        pychopper: sample_record[0].preprocessing == 'pychopper'
        bypass: sample_record[0].preprocessing == 'none'
    }

    PYCHOPPER(preprocessing_routes.pychopper.map { meta, reads, _primers -> tuple(meta, reads) })

    pychopper_reads = PYCHOPPER.out.processed_reads.map { meta, full_length, rescued ->
        tuple(meta, [full_length, rescued])
    }
    bypass_reads = preprocessing_routes.bypass.map { meta, reads, _primers ->
        tuple(meta, [reads])
    }
    alignment_reads = pychopper_reads.mix(bypass_reads)

    index_trigger = samples
        .map { _meta, _reads, _primers -> true }
        .first()
        .combine(fasta)
        .map { _trigger, reference -> reference }
    MINIMAP2_INDEX(index_trigger)
    ont_index = MINIMAP2_INDEX.out.index.first()

    MINIMAP2_SPLICE(alignment_reads, ont_index)

    versions = NANOPLOT_RAW_ONT.out.versions
        .mix(PYCHOPPER.out.versions)
        .mix(MINIMAP2_SPLICE.out.versions)

    emit:
    alignment = MINIMAP2_SPLICE.out.alignment
    raw_qc = NANOPLOT_RAW_ONT.out.qc
    preprocessing_reports = PYCHOPPER.out.reports
    versions = versions
}
