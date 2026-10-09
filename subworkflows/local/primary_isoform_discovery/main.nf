include { BLOCK2_CONTRACT_VALIDATE } from '../../../modules/local/block2_contract_validate/main'
include { ISOQUANT }                from '../../../modules/local/isoquant/main'
include { SQANTI3_QC }              from '../../../modules/local/sqanti3_qc/main'
include { TRANSCRIPTOME_CONTRACT }  from '../../../modules/local/transcriptome_contract/main'

workflow PRIMARY_ISOFORM_DISCOVERY {
    take:
    normalized_alignments
    reference_validation
    fasta
    gtf
    genome_build
    annotation_release

    main:
    BLOCK2_CONTRACT_VALIDATE(
        normalized_alignments,
        reference_validation,
        fasta,
        gtf,
        genome_build,
        annotation_release
    )

    ISOQUANT(BLOCK2_CONTRACT_VALIDATE.out.validated, fasta, gtf)
    SQANTI3_QC(ISOQUANT.out.results, fasta, gtf)

    keyed_isoquant_versions = ISOQUANT.out.versions.map { meta, versions ->
        tuple(meta.id, meta, versions)
    }
    keyed_sqanti3_versions = SQANTI3_QC.out.versions.map { meta, versions ->
        tuple(meta.id, versions)
    }
    paired_versions = keyed_isoquant_versions
        .join(keyed_sqanti3_versions)
        .map { sample_id, meta, isoquant_versions, sqanti3_versions ->
            tuple(sample_id, isoquant_versions, sqanti3_versions)
        }

    keyed_sqanti_results = SQANTI3_QC.out.results.map {
        meta, isoquant_dir, primary_gtf, support, counts, sqanti_dir, classification,
        junctions, transcript_fasta, block1_manifest, input_validation ->
        tuple(meta.id, meta, isoquant_dir, primary_gtf, support, counts, sqanti_dir,
              classification, junctions, transcript_fasta, block1_manifest, input_validation)
    }
    contract_inputs = keyed_sqanti_results
        .join(paired_versions)
        .map {
            _sample_id, meta, isoquant_dir, primary_gtf, support, counts, sqanti_dir,
            classification, junctions, transcript_fasta, block1_manifest, input_validation,
            isoquant_versions, sqanti3_versions ->
            tuple(meta, isoquant_dir, primary_gtf, support, counts, sqanti_dir,
                  classification, junctions, transcript_fasta, block1_manifest, input_validation,
                  isoquant_versions, sqanti3_versions)
        }

    TRANSCRIPTOME_CONTRACT(
        contract_inputs,
        reference_validation,
        genome_build,
        annotation_release
    )

    emit:
    transcriptome = TRANSCRIPTOME_CONTRACT.out.transcriptome
    summaries = TRANSCRIPTOME_CONTRACT.out.human_summary
    isoquant = ISOQUANT.out.results
    sqanti3 = SQANTI3_QC.out.results
    versions = ISOQUANT.out.versions.mix(SQANTI3_QC.out.versions)
}
