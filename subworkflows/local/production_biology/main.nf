include { SAMPLE_ABUNDANCE }       from '../../../modules/local/sample_abundance/main'
include { ABUNDANCE_MATRIX }        from '../../../modules/local/abundance_matrix/main'
include { EXPRESSION_QC }           from '../../../modules/local/expression_qc/main'
include { EDGER_DIFFERENTIAL }      from '../../../modules/local/edger_differential/main'
include { SATURN_DTU }              from '../../../modules/local/saturn_dtu/main'
include { SUPPA_ANALYSIS }          from '../../../modules/local/suppa_analysis/main'
include { TRANSCRIPT_ENDS }         from '../../../modules/local/transcript_ends/main'
include { TRANSDECODER }            from '../../../modules/local/transdecoder/main'
include { ISOFORM_SWITCH_ANALYZE }  from '../../../modules/local/isoform_switch_analyze/main'
include { FUNCTIONAL_ENRICHMENT }   from '../../../modules/local/functional_enrichment/main'
include { BIOLOGY_INTEGRATE }       from '../../../modules/local/biology_integrate/main'
include { BIOLOGY_CONTRACT }        from '../../../modules/local/biology_contract/main'

workflow PRODUCTION_BIOLOGY {
    take:
    transcriptomes
    isoquant_results

    main:
    keyed_transcriptomes = transcriptomes.map { meta, gtf, fasta, classification, support, sqanti, summary, structure, contract ->
        tuple(meta.id, meta, gtf, fasta, classification, support, sqanti, summary, structure, contract)
    }
    keyed_isoquant = isoquant_results.map { meta, isoquant_dir, _primary_gtf, _support, _counts, _manifest, _validation ->
        tuple(meta.id, isoquant_dir)
    }
    sample_inputs = keyed_transcriptomes.join(keyed_isoquant).map { _id, meta, gtf, fasta, classification, support, sqanti, summary, structure, contract, isoquant_dir ->
        tuple(meta, gtf, fasta, classification, support, sqanti, summary, structure, contract, isoquant_dir)
    }
    SAMPLE_ABUNDANCE(sample_inputs)
    bundle_dirs = SAMPLE_ABUNDANCE.out.bundle.map { _meta, bundle -> bundle }.collect()
    ABUNDANCE_MATRIX(bundle_dirs, params.contrast_numerator.toString(), params.contrast_denominator.toString())
    EXPRESSION_QC(ABUNDANCE_MATRIX.out.dataset)
    EDGER_DIFFERENTIAL(ABUNDANCE_MATRIX.out.dataset, params.contrast_numerator.toString(), params.contrast_denominator.toString(), params.block4_design_formula.toString())
    SATURN_DTU(ABUNDANCE_MATRIX.out.dataset, params.contrast_numerator.toString(), params.contrast_denominator.toString(), params.block4_design_formula.toString())
    SUPPA_ANALYSIS(ABUNDANCE_MATRIX.out.dataset, params.contrast_numerator.toString(), params.contrast_denominator.toString())
    TRANSCRIPT_ENDS(ABUNDANCE_MATRIX.out.dataset)
    TRANSDECODER(ABUNDANCE_MATRIX.out.dataset)
    ISOFORM_SWITCH_ANALYZE(ABUNDANCE_MATRIX.out.dataset, params.contrast_numerator.toString(), params.contrast_denominator.toString())
    FUNCTIONAL_ENRICHMENT(ABUNDANCE_MATRIX.out.dataset, EDGER_DIFFERENTIAL.out.results, SATURN_DTU.out.results, ISOFORM_SWITCH_ANALYZE.out.results)
    BIOLOGY_INTEGRATE(ABUNDANCE_MATRIX.out.dataset, EDGER_DIFFERENTIAL.out.results, SATURN_DTU.out.results, SUPPA_ANALYSIS.out.results, TRANSCRIPT_ENDS.out.results, TRANSDECODER.out.results, ISOFORM_SWITCH_ANALYZE.out.results, FUNCTIONAL_ENRICHMENT.out.results)
    BIOLOGY_CONTRACT(ABUNDANCE_MATRIX.out.dataset, EXPRESSION_QC.out.results, EDGER_DIFFERENTIAL.out.results, SATURN_DTU.out.results, SUPPA_ANALYSIS.out.results, TRANSCRIPT_ENDS.out.results, TRANSDECODER.out.results, ISOFORM_SWITCH_ANALYZE.out.results, FUNCTIONAL_ENRICHMENT.out.results, BIOLOGY_INTEGRATE.out.results)

    emit:
    contract = BIOLOGY_CONTRACT.out.contract
    dataset = ABUNDANCE_MATRIX.out.dataset
    integrated = BIOLOGY_INTEGRATE.out.results
    versions = EDGER_DIFFERENTIAL.out.versions.mix(SATURN_DTU.out.versions, SUPPA_ANALYSIS.out.versions, TRANSDECODER.out.versions, ISOFORM_SWITCH_ANALYZE.out.versions, FUNCTIONAL_ENRICHMENT.out.versions)
}
