include { CALLER_BENCHMARK } from '../subworkflows/local/caller_benchmark/main'

workflow BLOCK3_BENCHMARK {
    take:
    normalized_alignments
    reference_validation
    fasta
    gtf
    genome_build
    annotation_release

    main:
    CALLER_BENCHMARK(normalized_alignments, reference_validation, fasta, gtf, genome_build, annotation_release)

    emit:
    benchmark = CALLER_BENCHMARK.out.benchmark
    isoquant = CALLER_BENCHMARK.out.isoquant
    flair = CALLER_BENCHMARK.out.flair
    bambu = CALLER_BENCHMARK.out.bambu
    versions = CALLER_BENCHMARK.out.versions
}
