#!/usr/bin/env nextflow

nextflow.enable.dsl = 2

include { BLOCK1_ALIGNMENT } from './workflows/block1_alignment'
include { BLOCK2_TRANSCRIPTOME } from './workflows/block2_transcriptome'
include { BLOCK3_BENCHMARK } from './workflows/block3_benchmark'
include { BLOCK4_BIOLOGY } from './workflows/block4_biology'
include { FINAL_REPORTING } from './subworkflows/local/final_reporting/main'

def parseBooleanField(value, String field, String sample, boolean defaultValue = false) {
    def normalized = value?.toString()?.trim()?.toLowerCase()
    if (!normalized) {
        return defaultValue
    }
    if (!(normalized in ['true', 'false'])) {
        error "Sample '${sample}' has invalid ${field} value '${value}'. Expected 'true' or 'false'."
    }
    normalized == 'true'
}

def validateLabel(value, String field, String sample) {
    if (value.contains('\t') || value.contains('\n') || value.contains('\r')) {
        error "Sample '${sample}' has invalid ${field}: tabs and line breaks are not allowed."
    }
    value
}

workflow {
    if (params.help) {
        log.info """
        longread-rnaseq-pipeline ${workflow.manifest.version}

        Usage:
          nextflow run . -profile docker --input samples.csv --fasta genome.fa --gtf annotation.gtf \\
            --genome_build GRCh38 --annotation_release 'GENCODE release' [options]

        Required:
          --input                 CSV sample sheet (see docs/input.md)
          --fasta                 Reference genome FASTA
          --gtf                   Matching reference annotation GTF
          --genome_build          Declared genome build
          --annotation_release    Declared annotation release

        Modes:
          --mode default          Production IsoQuant + SQANTI3 (default)
          --mode benchmark        Default reconstruction plus FLAIR/Bambu comparison
          --mode biology          Production IsoQuant downstream biological analysis

        Common options:
          --outdir DIR            Output directory (default: results)
          --report_run_id ID      Portable run identifier
          --contrast_numerator C  Biology-mode contrast numerator
          --contrast_denominator C Biology-mode contrast denominator
          --help                  Show this message

        Configuration errors are checked before processes start. PacBio Iso-Seq
        HiFi rows require primer_fasta. Differential testing requires an explicit
        contrast and adequate biological replication. Full parameter definitions
        and supported sample-sheet combinations are in nextflow_schema.json and
        docs/input.md.
        """.stripIndent().trim()
        return
    }

    def required_parameters = [
        'input',
        'fasta',
        'gtf',
        'genome_build',
        'annotation_release'
    ]
    def missing_parameters = required_parameters.findAll { key ->
        !params[key]?.toString()?.trim()
    }
    if (missing_parameters) {
        error "Missing required parameter(s): ${missing_parameters.collect { parameter -> '--' + parameter }.join(', ')}"
    }

    if (!(params.mode in ['default', 'benchmark', 'biology'])) {
        error "Invalid --mode '${params.mode}'. Expected 'default', 'benchmark', or 'biology'."
    }
    if (params.annotation_guided_alignment) {
        error "--annotation_guided_alignment is reserved but not implemented in Block 1B. The default unguided alignment path must remain active."
    }
    if (!(params.large_file_publish_mode in ['copy', 'link', 'rellink', 'symlink'])) {
        error "Invalid --large_file_publish_mode '${params.large_file_publish_mode}'. Expected copy, link, rellink, or symlink."
    }
    if (params.transcriptome_filtering) {
        error "--transcriptome_filtering is reserved for an explicit future curation policy and is not implemented in Block 2. Raw IsoQuant models remain unfiltered."
    }
    if (params.isoquant_report_novel_unspliced != null && !(params.isoquant_report_novel_unspliced in [true, false])) {
        error "--isoquant_report_novel_unspliced must be true, false, or omitted to use IsoQuant's platform-specific default."
    }
    if (params.flair_min_support < 1 || params.flair_end_window < 0) {
        error "--flair_min_support must be at least 1 and --flair_end_window must be non-negative."
    }
    if (params.benchmark_single_exon_reciprocal_overlap <= 0 || params.benchmark_single_exon_reciprocal_overlap > 1) {
        error "--benchmark_single_exon_reciprocal_overlap must be greater than 0 and at most 1."
    }
    if ((params.contrast_numerator && !params.contrast_denominator) || (!params.contrast_numerator && params.contrast_denominator)) {
        error "--contrast_numerator and --contrast_denominator must be supplied together."
    }
    if (params.contrast_numerator && params.contrast_numerator == params.contrast_denominator) {
        error "Block 4 contrast numerator and denominator must be distinct."
    }
    if (params.mode != 'biology' && (params.contrast_numerator || params.contrast_denominator)) {
        error "Differential contrasts are valid only with --mode biology. Remove the contrast or select biology mode."
    }
    if (params.mode != 'biology' && params.enable_functional_enrichment) {
        error "--enable_functional_enrichment requires --mode biology."
    }
    if (params.report_run_id && !(params.report_run_id.toString() ==~ /[A-Za-z0-9_.-]+/)) {
        error "--report_run_id may contain only letters, numbers, '.', '_', and '-'."
    }
    if (params.max_cpus < 1) {
        error "--max_cpus must be at least 1."
    }
    if (params.minimum_biological_replicates < 2) {
        error "--minimum_biological_replicates must be at least 2."
    }
    if (!params.block4_design_formula.toString().contains('condition') || !(params.block4_design_formula.toString() =~ /(^|[^0-9])0([^0-9]|$)/)) {
        error "--block4_design_formula must include condition and use a no-intercept (0+) parameterization."
    }
    if (params.block4_fdr <= 0 || params.block4_fdr > 1 || params.isoform_switch_min_dif < 0 || params.isoform_switch_min_dif > 1) {
        error "--block4_fdr must be in (0,1] and --isoform_switch_min_dif must be in [0,1]."
    }
    if (params.transcript_end_cluster_window < 0 || params.transdecoder_min_peptide_length < 1) {
        error "Block 4 coordinate windows must be non-negative and minimum peptide length must be positive."
    }
    def benchmark_end_bins = params.benchmark_end_bins.toString().split(',').collect { value -> value as Integer }
    if (!benchmark_end_bins || benchmark_end_bins[0] != 0 || benchmark_end_bins != benchmark_end_bins.sort().unique()) {
        error "--benchmark_end_bins must be an ascending, unique comma-separated list beginning with 0."
    }

    def reference_fasta = file(params.fasta, checkIfExists: true)
    def annotation_gtf = file(params.gtf, checkIfExists: true)

    log.info """
    longread-rnaseq-pipeline ${workflow.manifest.version} (Blocks 1-5)
    ----------------------------------------------------------------
    Sample sheet                : ${params.input}
    Output directory            : ${params.outdir}
    Genome build                : ${params.genome_build}
    Reference FASTA             : ${reference_fasta}
    Annotation                  : ${annotation_gtf}
    Annotation release          : ${params.annotation_release}
    Pipeline mode               : ${params.mode}
    Annotation-guided alignment : ${params.annotation_guided_alignment}
    PacBio cluster2 singletons  : ${params.isoseq_include_singletons}
    IsoQuant novel unspliced    : ${params.isoquant_report_novel_unspliced == null ? 'tool-native by data type' : params.isoquant_report_novel_unspliced}
    Transcriptome filtering     : ${params.transcriptome_filtering}
    Execution profile           : ${params.execution_profile}
    ----------------------------------------------------------------
    Block 1 performs platform-specific preprocessing, alignment, and QC.
    Default mode continues with IsoQuant discovery and SQANTI3 structural QC.
    Benchmark mode compares IsoQuant, FLAIR, and Bambu through common SQANTI3 characterization.
    Biology mode uses only the production IsoQuant transcriptome for downstream biological analyses.
    Every mode ends with deterministic QC aggregation, final reporting, and consolidated provenance.
    """.stripIndent().trim()

    samplesheet = channel
        .fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, strip: true)
        .ifEmpty { error "The sample sheet contains no sample records: ${params.input}" }

    validated_samples = samplesheet.map { row ->
        def columns = [
            'sample',
            'platform',
            'library_type',
            'preprocessing',
            'condition',
            'replicate',
            'reads',
            'primer_fasta',
            'primer_pair',
            'ont_kit',
            'require_polya'
        ]
        def absent_columns = columns.findAll { column -> !row.containsKey(column) }
        if (absent_columns) {
            error "Sample sheet is missing required column(s): ${absent_columns.join(', ')}"
        }

        def record = columns.collectEntries { column ->
            [(column): row[column]?.toString()?.trim()]
        }
        def required_values = ['sample', 'platform', 'library_type', 'preprocessing', 'condition', 'replicate', 'reads']
        def empty_values = required_values.findAll { column -> !record[column] }
        if (empty_values) {
            error "Sample '${record.sample ?: '<unknown>'}' has empty required field(s): ${empty_values.join(', ')}"
        }
        if (!(record.sample ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
            error "Sample '${record.sample}' is invalid. Use letters, numbers, '.', '_', or '-', beginning with a letter or number."
        }

        validateLabel(record.condition, 'condition', record.sample)
        validateLabel(record.replicate, 'replicate', record.sample)

        def supported_library_types = [
            ont: ['cdna', 'pcr_cdna'],
            pacbio: ['isoseq_hifi', 'isoseq_flnc']
        ]
        if (!supported_library_types.containsKey(record.platform) ||
            !(record.library_type in supported_library_types[record.platform])) {
            error "Unsupported platform/library combination '${record.platform}/${record.library_type}' for sample '${record.sample}'. Supported combinations: ont/cdna, ont/pcr_cdna, pacbio/isoseq_hifi, pacbio/isoseq_flnc."
        }

        def reads_lower = record.reads.toLowerCase()
        if (record.platform == 'ont' && !(reads_lower ==~ /.*\.(fastq|fq)(\.gz)?$/)) {
            error "ONT sample '${record.sample}' requires FASTQ input ending in .fastq, .fq, .fastq.gz, or .fq.gz."
        }
        if (record.platform == 'pacbio' && !reads_lower.endsWith('.bam')) {
            error "PacBio sample '${record.sample}' requires BAM input."
        }

        def require_polya = parseBooleanField(record.require_polya, 'require_polya', record.sample, false)
        def primer_path = []

        if (record.platform == 'ont') {
            if (!(record.preprocessing in ['pychopper', 'none'])) {
                error "ONT sample '${record.sample}' requires preprocessing 'pychopper' or 'none'."
            }
            if (record.preprocessing == 'pychopper') {
                if (!(record.ont_kit in ['PCS109', 'PCS110', 'PCS111', 'PCS114', 'PCB111', 'PCB114'])) {
                    error "ONT sample '${record.sample}' uses Pychopper and requires ont_kit PCS109, PCS110, PCS111, PCS114, PCB111, or PCB114."
                }
            } else if (record.ont_kit) {
                error "ONT sample '${record.sample}' sets ont_kit but preprocessing is 'none'."
            }
            if (record.primer_fasta) {
                error "ONT sample '${record.sample}' must not set primer_fasta; that field is for PacBio Iso-Seq processing."
            }
            if (record.primer_pair) {
                error "ONT sample '${record.sample}' must not set primer_pair; that field selects a PacBio lima output."
            }
            if (require_polya) {
                error "ONT sample '${record.sample}' must not enable require_polya; it controls PacBio isoseq refine."
            }
        }

        if (record.platform == 'pacbio' && record.library_type == 'isoseq_hifi') {
            if (record.preprocessing != 'isoseq') {
                error "PacBio HiFi sample '${record.sample}' requires preprocessing 'isoseq' (lima, refine, cluster2)."
            }
            if (!record.primer_fasta) {
                error "PacBio HiFi sample '${record.sample}' requires an explicit primer_fasta."
            }
            primer_path = file(record.primer_fasta, checkIfExists: true)
            if (record.primer_pair && !(record.primer_pair ==~ /[A-Za-z0-9][A-Za-z0-9_.-]*--[A-Za-z0-9][A-Za-z0-9_.-]*/)) {
                error "PacBio HiFi sample '${record.sample}' has invalid primer_pair '${record.primer_pair}'. Expected a lima pair such as IsoSeqX_bc06_5p--IsoSeqX_3p."
            }
            if (record.ont_kit) {
                error "PacBio sample '${record.sample}' must not set ont_kit."
            }
        }

        if (record.platform == 'pacbio' && record.library_type == 'isoseq_flnc') {
            if (record.preprocessing != 'none') {
                error "PacBio FLNC sample '${record.sample}' requires preprocessing 'none'; lima/refine are bypassed and cluster2 remains active."
            }
            if (record.primer_fasta) {
                error "PacBio FLNC sample '${record.sample}' must not set primer_fasta because lima/refine are bypassed."
            }
            if (record.primer_pair) {
                error "PacBio FLNC sample '${record.sample}' must not set primer_pair because lima is bypassed."
            }
            if (record.ont_kit) {
                error "PacBio sample '${record.sample}' must not set ont_kit."
            }
            if (require_polya) {
                error "PacBio FLNC sample '${record.sample}' cannot enable require_polya because refine is bypassed."
            }
        }

        def meta = [
            id: record.sample,
            platform: record.platform,
            library_type: record.library_type,
            preprocessing: record.preprocessing,
            condition: record.condition,
            replicate: record.replicate,
            original_input: record.reads,
            primer_fasta: record.primer_fasta ?: null,
            primer_pair: record.primer_pair ?: null,
            ont_kit: record.ont_kit ?: null,
            require_polya: require_polya
        ]
        tuple(meta, file(record.reads, checkIfExists: true), primer_path)
    }

    checked_samples = validated_samples
        .collect(flat: false)
        .flatMap { records ->
            def duplicate_ids = records
                .countBy { item -> item[0].id }
                .findAll { _sample_id, count -> count > 1 }
                .keySet()
            if (duplicate_ids) {
                error "Sample identifiers must be unique. Duplicates: ${duplicate_ids.sort().join(', ')}"
            }
            records
        }

    BLOCK1_ALIGNMENT(
        checked_samples,
        channel.value(reference_fasta),
        channel.value(reference_fasta.toString()),
        channel.value(annotation_gtf.toString()),
        channel.value(params.genome_build.toString()),
        channel.value(params.annotation_release.toString())
    )

    BLOCK1_ALIGNMENT.out.normalized.view { meta, _bam, _bai, _qc_dir, _manifest ->
        "Normalized aligned-read record ready: ${meta.id} (${meta.platform}/${meta.library_type})"
    }

    def report_placeholder = file("${projectDir}/assets/reporting/not_run.json", checkIfExists: true)
    normalized_contracts = BLOCK1_ALIGNMENT.out.normalized.map { _meta, _bam, _bai, _qc, manifest -> manifest }
    qc_inputs = BLOCK1_ALIGNMENT.out.alignment_qc.map { _meta, qc -> qc }
        .mix(BLOCK1_ALIGNMENT.out.ont_raw_qc.map { _meta, qc -> qc })
        .mix(BLOCK1_ALIGNMENT.out.pacbio_raw_qc.map { _meta, qc -> qc })
    transcriptome_contracts = channel.of(report_placeholder)
    benchmark_contracts = channel.of(report_placeholder)
    biology_contracts = channel.of(report_placeholder)
    transcript_summaries = channel.of(report_placeholder)
    analysis_versions = channel.empty()

    if (params.mode in ['default', 'biology']) {
        BLOCK2_TRANSCRIPTOME(
            BLOCK1_ALIGNMENT.out.normalized,
            BLOCK1_ALIGNMENT.out.reference_validation,
            channel.value(reference_fasta),
            channel.value(annotation_gtf),
            channel.value(params.genome_build.toString()),
            channel.value(params.annotation_release.toString())
        )

        BLOCK2_TRANSCRIPTOME.out.transcriptome.view { meta, _gtf, _fasta, _classification, _support, _sqanti, _summary, _structure, _contract ->
            "Standardized transcriptome record ready: ${meta.id} (${meta.platform}/${meta.library_type})"
        }
        transcriptome_contracts = BLOCK2_TRANSCRIPTOME.out.transcriptome.map { _meta, _gtf, _fasta, _classification, _support, _sqanti, summary, _structure, contract -> contract }
        transcript_summaries = BLOCK2_TRANSCRIPTOME.out.transcriptome.map { _meta, _gtf, _fasta, _classification, _support, _sqanti, summary, _structure, _contract -> summary }
        analysis_versions = BLOCK2_TRANSCRIPTOME.out.versions.map { _meta, version -> version }
        if (params.mode == 'biology') {
            BLOCK4_BIOLOGY(BLOCK2_TRANSCRIPTOME.out.transcriptome, BLOCK2_TRANSCRIPTOME.out.isoquant)
            BLOCK4_BIOLOGY.out.contract.view { contract -> "Block 4 biology contract ready: ${contract}" }
            biology_contracts = BLOCK4_BIOLOGY.out.contract
            analysis_versions = analysis_versions.mix(BLOCK4_BIOLOGY.out.versions)
        }
    } else {
        BLOCK3_BENCHMARK(
            BLOCK1_ALIGNMENT.out.normalized,
            BLOCK1_ALIGNMENT.out.reference_validation,
            channel.value(reference_fasta),
            channel.value(annotation_gtf),
            channel.value(params.genome_build.toString()),
            channel.value(params.annotation_release.toString())
        )

        BLOCK3_BENCHMARK.out.benchmark.view { meta, _summary, _contract, _membership, _pairwise, _ends, _partial, _genes ->
            "Caller benchmark record ready: ${meta.id} (${meta.platform}/${meta.library_type})"
        }
        benchmark_contracts = BLOCK3_BENCHMARK.out.benchmark.map { _meta, _summary, contract, _membership, _pairwise, _ends, _partial, _genes -> contract }
        analysis_versions = BLOCK3_BENCHMARK.out.versions.map { _meta, version -> version }
    }

    def resolved_parameters = [
        input: params.input, outdir: params.outdir, fasta: params.fasta, gtf: params.gtf,
        genome_build: params.genome_build, annotation_release: params.annotation_release,
        mode: params.mode, execution_profile: params.execution_profile,
        annotation_guided_alignment: params.annotation_guided_alignment,
        lima_peek_guess: params.lima_peek_guess,
        isoseq_include_singletons: params.isoseq_include_singletons,
        large_file_publish_mode: params.large_file_publish_mode,
        isoquant_report_novel_unspliced: params.isoquant_report_novel_unspliced,
        transcriptome_filtering: params.transcriptome_filtering,
        flair_min_support: params.flair_min_support, flair_end_window: params.flair_end_window,
        flair_parallel_mode: params.flair_parallel_mode, benchmark_end_bins: params.benchmark_end_bins,
        benchmark_single_exon_reciprocal_overlap: params.benchmark_single_exon_reciprocal_overlap,
        contrast_numerator: params.contrast_numerator, contrast_denominator: params.contrast_denominator,
        block4_design_formula: params.block4_design_formula,
        minimum_biological_replicates: params.minimum_biological_replicates,
        block4_fdr: params.block4_fdr, transcript_end_cluster_window: params.transcript_end_cluster_window,
        transdecoder_min_peptide_length: params.transdecoder_min_peptide_length,
        transdecoder_single_best_only: params.transdecoder_single_best_only,
        isoform_switch_min_dif: params.isoform_switch_min_dif,
        enable_functional_enrichment: params.enable_functional_enrichment,
        max_cpus: params.max_cpus, max_memory: params.max_memory, max_time: params.max_time
    ]
    all_versions = BLOCK1_ALIGNMENT.out.runtime_versions.map { _meta, version -> version }
        .mix(BLOCK1_ALIGNMENT.out.reference_versions)
        .mix(analysis_versions)
    FINAL_REPORTING(
        normalized_contracts,
        transcriptome_contracts,
        benchmark_contracts,
        biology_contracts,
        qc_inputs,
        transcript_summaries,
        all_versions,
        BLOCK1_ALIGNMENT.out.reference_validation,
        channel.value(file(params.input, checkIfExists: true)),
        channel.value(file("${projectDir}/assets/reporting/status_vocabulary.json", checkIfExists: true)),
        channel.value(file("${projectDir}/assets/reporting/tool_registry.json", checkIfExists: true)),
        channel.value(file("${projectDir}/assets/reporting/multiqc_config.yaml", checkIfExists: true)),
        groovy.json.JsonOutput.toJson(resolved_parameters),
        params.mode.toString(),
        workflow.manifest.version.toString(),
        workflow.nextflow.version.toString(),
        (params.report_run_id?.toString()?.trim() ?: workflow.runName.toString()),
        params.execution_profile.toString(),
        workflow.start.toString()
    )
    FINAL_REPORTING.out.report.view { report -> "Final Block 5 report ready: ${report}" }

}
