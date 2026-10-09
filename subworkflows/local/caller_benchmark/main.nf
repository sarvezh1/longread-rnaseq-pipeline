include { BLOCK2_CONTRACT_VALIDATE }           from '../../../modules/local/block2_contract_validate/main'
include { ISOQUANT }                           from '../../../modules/local/isoquant/main'
include { FLAIR_TRANSCRIPTOME }                from '../../../modules/local/flair_transcriptome/main'
include { BAMBU }                              from '../../../modules/local/bambu/main'
include { SQANTI3_QC as SQANTI3_ISOQUANT }     from '../../../modules/local/sqanti3_qc/main'
include { SQANTI3_QC as SQANTI3_FLAIR }        from '../../../modules/local/sqanti3_qc/main'
include { SQANTI3_QC as SQANTI3_BAMBU }        from '../../../modules/local/sqanti3_qc/main'
include { TRANSCRIPT_STANDARDIZE as STANDARDIZE_ISOQUANT } from '../../../modules/local/transcript_standardize/main'
include { TRANSCRIPT_STANDARDIZE as STANDARDIZE_FLAIR }    from '../../../modules/local/transcript_standardize/main'
include { TRANSCRIPT_STANDARDIZE as STANDARDIZE_BAMBU }    from '../../../modules/local/transcript_standardize/main'
include { TRANSCRIPT_CONCORDANCE }             from '../../../modules/local/transcript_concordance/main'

workflow CALLER_BENCHMARK {
    take:
    normalized_alignments
    reference_validation
    fasta
    gtf
    genome_build
    annotation_release

    main:
    benchmark_inputs = normalized_alignments.map { meta, bam, bai, qc, manifest ->
        tuple(meta + [benchmark_mode: true], bam, bai, qc, manifest)
    }
    BLOCK2_CONTRACT_VALIDATE(benchmark_inputs, reference_validation, fasta, gtf, genome_build, annotation_release)

    ISOQUANT(BLOCK2_CONTRACT_VALIDATE.out.validated, fasta, gtf)
    FLAIR_TRANSCRIPTOME(BLOCK2_CONTRACT_VALIDATE.out.validated, fasta)
    BAMBU(BLOCK2_CONTRACT_VALIDATE.out.validated, fasta, gtf)

    isoquant_sqanti_input = ISOQUANT.out.results.map { meta, caller_dir, model_gtf, support, counts, manifest, validation ->
        def benchmark_meta = meta + [id: "${meta.id}.isoquant", benchmark_sample_id: meta.id, benchmark_caller: 'isoquant']
        tuple(benchmark_meta, caller_dir, model_gtf, support, counts, manifest, validation)
    }
    flair_sqanti_input = FLAIR_TRANSCRIPTOME.out.results.map { meta, caller_dir, model_gtf, support, counts, manifest, validation ->
        def benchmark_meta = meta + [id: "${meta.id}.flair", benchmark_sample_id: meta.id, benchmark_caller: 'flair']
        tuple(benchmark_meta, caller_dir, model_gtf, support, counts, manifest, validation)
    }
    bambu_sqanti_input = BAMBU.out.results.map { meta, caller_dir, model_gtf, support, counts, manifest, validation ->
        def benchmark_meta = meta + [id: "${meta.id}.bambu", benchmark_sample_id: meta.id, benchmark_caller: 'bambu']
        tuple(benchmark_meta, caller_dir, model_gtf, support, counts, manifest, validation)
    }

    SQANTI3_ISOQUANT(isoquant_sqanti_input, fasta, gtf)
    SQANTI3_FLAIR(flair_sqanti_input, fasta, gtf)
    SQANTI3_BAMBU(bambu_sqanti_input, fasta, gtf)

    STANDARDIZE_ISOQUANT(SQANTI3_ISOQUANT.out.results)
    STANDARDIZE_FLAIR(SQANTI3_FLAIR.out.results)
    STANDARDIZE_BAMBU(SQANTI3_BAMBU.out.results)

    keyed_isoquant = STANDARDIZE_ISOQUANT.out.standardized.map { meta, structures, junctions, summary, manifest ->
        tuple(meta.benchmark_sample_id, meta, structures, manifest)
    }
    keyed_flair = STANDARDIZE_FLAIR.out.standardized.map { meta, structures, junctions, summary, manifest ->
        tuple(meta.benchmark_sample_id, structures)
    }
    keyed_bambu = STANDARDIZE_BAMBU.out.standardized.map { meta, structures, junctions, summary, manifest ->
        tuple(meta.benchmark_sample_id, structures)
    }
    structures = keyed_isoquant.join(keyed_flair).join(keyed_bambu)

    iq_versions = ISOQUANT.out.versions.map { meta, version -> tuple(meta.id, version) }
    fl_versions = FLAIR_TRANSCRIPTOME.out.versions.map { meta, version -> tuple(meta.id, version) }
    ba_versions = BAMBU.out.versions.map { meta, version -> tuple(meta.id, version) }
    sq_iq_versions = SQANTI3_ISOQUANT.out.versions.map { meta, version -> tuple(meta.benchmark_sample_id, version) }
    sq_fl_versions = SQANTI3_FLAIR.out.versions.map { meta, version -> tuple(meta.benchmark_sample_id, version) }
    sq_ba_versions = SQANTI3_BAMBU.out.versions.map { meta, version -> tuple(meta.benchmark_sample_id, version) }
    version_bundle = iq_versions.join(fl_versions).join(ba_versions).join(sq_iq_versions).join(sq_fl_versions).join(sq_ba_versions)

    concordance_input = structures.join(version_bundle).map {
        sample_id, caller_meta, isoquant_structures, manifest, flair_structures, bambu_structures,
        isoquant_version, flair_version, bambu_version, sqanti_isoquant_version, sqanti_flair_version, sqanti_bambu_version ->
        def sample_meta = caller_meta + [id: sample_id]
        tuple(sample_meta, isoquant_structures, flair_structures, bambu_structures, manifest,
              isoquant_version, flair_version, bambu_version,
              sqanti_isoquant_version, sqanti_flair_version, sqanti_bambu_version)
    }
    TRANSCRIPT_CONCORDANCE(concordance_input)

    emit:
    benchmark = TRANSCRIPT_CONCORDANCE.out.benchmark
    isoquant = ISOQUANT.out.results
    flair = FLAIR_TRANSCRIPTOME.out.results
    bambu = BAMBU.out.results
    versions = TRANSCRIPT_CONCORDANCE.out.versions
}
