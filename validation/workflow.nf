nextflow.enable.dsl = 2

params.input = null
params.outdir = 'validation/results'
params.end_tolerance = 0

process VALIDATION_TRUTH_SCORE {
    tag "${meta.run_id}"
    cpus 1
    memory '2 GB'
    time '2 h'
    container 'quay.io/biocontainers/isoquant:4.0.0--pyh106432d_0@sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142'
    publishDir "${params.outdir}/sirv", mode: 'copy'

    input:
    tuple val(meta), path(prediction), path(truth)

    output:
    tuple val(meta), path("${meta.run_id}.truth_metrics.json"), path("${meta.run_id}.truth_matches.tsv"), emit: metrics

    script:
    def prefix = meta.truth_contig_prefix ? "--truth-contig-prefix '${meta.truth_contig_prefix}'" : ''
    """
    set -euo pipefail
    python ${projectDir}/../bin/validation_metrics.py truth-match \
      --prediction ${prediction} --prediction-format '${meta.prediction_format}' \
      --truth ${truth} ${prefix} --end-tolerance ${params.end_tolerance} \
      --dataset-id '${meta.dataset_id}' --caller '${meta.caller}' --platform '${meta.platform}' \
      --truth-type '${meta.truth_type}' --output ${meta.run_id}.truth_metrics.json \
      --matches ${meta.run_id}.truth_matches.tsv
    """

    stub:
    """
    printf '{"schema_version":"6.0","dataset_id":"${meta.dataset_id}","caller":"${meta.caller}","status":"STUB"}\n' > ${meta.run_id}.truth_metrics.json
    printf 'predicted_transcript_id\ttruth_transcript_id\n' > ${meta.run_id}.truth_matches.tsv
    """
}

workflow {
    if (!params.input) error 'Block 6 scoring requires --input validation scoring manifest.'
    if (params.end_tolerance < 0) error '--end_tolerance must be non-negative.'
    scoring = channel.fromPath(params.input, checkIfExists: true)
        .splitCsv(header: true, sep: '\t')
        .map { row ->
            def missing = ['run_id','dataset_id','platform','caller','truth_type','prediction','prediction_format','truth'].findAll { field -> !row[field] }
            if (missing) error "Validation scoring manifest fields are empty: ${missing.join(', ')}."
            if (!(row.prediction_format in ['gtf','standardized'])) error "Unsupported prediction_format '${row.prediction_format}'."
            def meta = [run_id: row.run_id, dataset_id: row.dataset_id, platform: row.platform,
                        caller: row.caller, truth_type: row.truth_type,
                        prediction_format: row.prediction_format,
                        truth_contig_prefix: row.truth_contig_prefix ?: '']
            tuple(meta, file(row.prediction, checkIfExists: true), file(row.truth, checkIfExists: true))
        }
    VALIDATION_TRUTH_SCORE(scoring)
}
