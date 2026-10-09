process LIMA_ISOSEQ {
    tag "${meta.id} (Iso-Seq primers)"
    label 'process_high'

    container 'quay.io/biocontainers/lima:26.2.1--h9ee0642_0'

    publishDir "${params.outdir}/preprocessing/pacbio",
        mode: params.large_file_publish_mode,
        pattern: "*.lima.bam*"
    publishDir "${params.outdir}/preprocessing/pacbio",
        mode: 'copy',
        pattern: "*.lima*.{log,report,counts,json,csv}"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.lima.versions.yml"

    input:
    tuple val(meta), path(reads), path(primers)

    output:
    tuple val(meta), path("${meta.id}.lima.bam"), path(primers), emit: full_length_bam
    tuple val(meta), path("${meta.id}.lima.log"), emit: log
    tuple val(meta), path("${meta.id}.lima.bam.pbi"), optional: true, emit: pbi
    tuple val(meta), path("*.lima.report"), optional: true, emit: report
    tuple val(meta), path("*.lima.counts"), optional: true, emit: counts
    tuple val(meta), path("${meta.id}.lima.versions.yml"), emit: versions

    script:
    def peek_guess = params.lima_peek_guess ? '--peek-guess' : ''
    """
    lima \
        --isoseq \
        ${peek_guess} \
        -j ${task.cpus} \
        ${reads} \
        ${primers} \
        ${meta.id}.lima_output.bam \
        > ${meta.id}.lima.log 2>&1

    mapfile -t lima_bams < <(find . -maxdepth 1 -type f -name '${meta.id}.lima_output*.bam' | sort)
    primer_pair='${meta.primer_pair ?: ''}'
    if [[ -n "\${primer_pair}" ]]; then
        source_bam='${meta.id}.lima_output.'"\${primer_pair}"'.bam'
        if [[ ! -s "\${source_bam}" ]]; then
            echo "Requested primer_pair '\${primer_pair}' did not match a non-empty lima output for sample ${meta.id}." >&2
            printf 'Available outputs:\\n' >&2
            printf '  %s\\n' "\${lima_bams[@]}" >&2
            exit 1
        fi
    elif [[ \${#lima_bams[@]} -eq 1 ]]; then
        source_bam="\${lima_bams[0]}"
    else
        echo "lima produced \${#lima_bams[@]} BAMs for sample ${meta.id}; set sample-sheet primer_pair to select exactly one output." >&2
        printf 'Available outputs:\\n' >&2
        printf '  %s\\n' "\${lima_bams[@]}" >&2
        exit 1
    fi

    mv "\${source_bam}" ${meta.id}.lima.bam
    if [[ -f "\${source_bam}.pbi" ]]; then
        mv "\${source_bam}.pbi" ${meta.id}.lima.bam.pbi
    fi

    cat > ${meta.id}.lima.versions.yml <<VERSIONS
    LIMA_ISOSEQ:
      lima: "\$(lima --version 2>&1 | sed -n '1s/^lima //p' | awk '{print \$1}')"
      isoseq_mode: true
      peek_guess: ${params.lima_peek_guess}
      selected_primer_pair: "${meta.primer_pair ?: 'automatic-single-output'}"
      container: "quay.io/biocontainers/lima:26.2.1--h9ee0642_0"
    VERSIONS
    """

    stub:
    """
    touch ${meta.id}.lima.bam
    touch ${meta.id}.lima.bam.pbi
    touch ${meta.id}.lima.log
    cat > ${meta.id}.lima.versions.yml <<VERSIONS
    LIMA_ISOSEQ:
      lima: "26.2.1"
      isoseq_mode: true
      peek_guess: ${params.lima_peek_guess}
      selected_primer_pair: "${meta.primer_pair ?: 'automatic-single-output'}"
      container: "quay.io/biocontainers/lima:26.2.1--h9ee0642_0"
    VERSIONS
    """
}
