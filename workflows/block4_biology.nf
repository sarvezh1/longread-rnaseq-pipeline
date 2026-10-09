include { PRODUCTION_BIOLOGY } from '../subworkflows/local/production_biology/main'

workflow BLOCK4_BIOLOGY {
    take:
    transcriptomes
    isoquant_results
    main:
    PRODUCTION_BIOLOGY(transcriptomes, isoquant_results)
    emit:
    contract = PRODUCTION_BIOLOGY.out.contract
    dataset = PRODUCTION_BIOLOGY.out.dataset
    integrated = PRODUCTION_BIOLOGY.out.integrated
    versions = PRODUCTION_BIOLOGY.out.versions
}
