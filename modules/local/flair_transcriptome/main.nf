process FLAIR_TRANSCRIPTOME {
    tag "${meta.id} (${meta.platform})"
    label 'process_high'

    container 'quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0@sha256:a39894c1507c68241fba62ddbe1a931d5a6bfaacbcb6be24631d46b60e67ce84'

    publishDir "${params.outdir}/benchmark/flair", mode: 'copy', pattern: '*.flair'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: '*.flair.versions.yml'

    input:
    tuple val(meta), path(bam), path(bai), path(alignment_qc), path(block1_manifest), path(input_validation)
    path fasta

    output:
    tuple val(meta), path("${meta.id}.flair"), path("${meta.id}.flair.gtf"), path("${meta.id}.flair.support.tsv"), path("${meta.id}.flair.counts.tsv"), path(block1_manifest), path(input_validation), emit: results
    tuple val(meta), path("${meta.id}.flair.versions.yml"), emit: versions

    script:
    def sample_id = meta.id
    def min_support = params.flair_min_support
    def end_window = params.flair_end_window
    def parallel_mode = params.flair_parallel_mode
    """
    set -euo pipefail
    mkdir -p ${sample_id}.flair
    set +e
    flair transcriptome \
        --genomealignedbam ${bam} \
        --genome ${fasta} \
        --output ${sample_id}.flair/${sample_id} \
        --threads ${task.cpus} \
        --support ${min_support} \
        --end_window ${end_window} \
        --parallelmode '${parallel_mode}' \
        --noaligntoannot \
        2>&1 | tee ${sample_id}.flair/execution.log
    flair_status=\${PIPESTATUS[0]}
    set -e

    gtf_file='${sample_id}.flair/${sample_id}.isoforms.gtf'
    bed_file='${sample_id}.flair/${sample_id}.isoforms.bed'
    fasta_file='${sample_id}.flair/${sample_id}.isoforms.fa'
    map_file='${sample_id}.flair/${sample_id}.isoform.read.map.txt'
    if [ "\${flair_status}" -ne 0 ]; then
        if grep -Fq 'No junctions from GTF or junctionsBed to correct with' ${sample_id}.flair/execution.log; then
            printf '# no FLAIR transcript models: normalized BAM had no usable splice junctions\n' > "\${gtf_file}"
            : > "\${bed_file}"
            : > "\${fasta_file}"
            : > "\${map_file}"
            execution_status='NO_USABLE_SPLICE_JUNCTIONS'
        else
            exit "\${flair_status}"
        fi
    else
        test -f "\${gtf_file}"
        test -f "\${bed_file}"
        test -f "\${fasta_file}"
        test -f "\${map_file}"
        execution_status='COMPLETE'
    fi
    cp "\${gtf_file}" ${sample_id}.flair.gtf
    awk -F '\t' 'BEGIN {OFS="\t"; print "transcript_id","count"} NF {n=split(\$2,a,","); print \$1,n}' "\${map_file}" > ${sample_id}.flair.support.tsv
    cp ${sample_id}.flair.support.tsv ${sample_id}.flair.counts.tsv
    cp ${sample_id}.flair.support.tsv ${sample_id}.flair/

    flair_version=\$(flair --version 2>&1 | awk 'NF {print \$NF; exit}')
    printf '%s\n' \
        'FLAIR:' \
        "  flair: \"\${flair_version}\"" \
        '  container: "quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0@sha256:a39894c1507c68241fba62ddbe1a931d5a6bfaacbcb6be24631d46b60e67ce84"' \
        '  workflow: "transcriptome"' \
        '  input: "coordinate-sorted indexed genome-aligned BAM"' \
        '  annotation_alignment: false' \
        '  reference_gtf_supplied: false' \
        '  min_support: ${min_support}' \
        '  end_window: ${end_window}' \
        '  parallel_mode: "${parallel_mode}"' \
        "  status: \"\${execution_status}\"" \
        > ${sample_id}.flair.versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}.flair
    printf 'chrStub\tFLAIR\ttranscript\t1\t100\t.\t+\t.\tgene_id "flair_gene"; transcript_id "flair_tx";\nchrStub\tFLAIR\texon\t1\t100\t.\t+\t.\tgene_id "flair_gene"; transcript_id "flair_tx";\n' > ${meta.id}.flair.gtf
    printf 'transcript_id\tcount\nflair_tx\t3\n' > ${meta.id}.flair.support.tsv
    cp ${meta.id}.flair.support.tsv ${meta.id}.flair.counts.tsv
    cp ${meta.id}.flair.gtf ${meta.id}.flair/
    printf '>flair_tx\nACGT\n' > ${meta.id}.flair/${meta.id}.isoforms.fasta
    printf 'flair_tx\tread1,read2,read3\n' > ${meta.id}.flair/${meta.id}.isoform.read.map.txt
    printf 'chrStub\t0\t100\tflair_tx\n' > ${meta.id}.flair/${meta.id}.isoforms.bed
    printf 'FLAIR:\n  flair: "3.0.1"\n' > ${meta.id}.flair.versions.yml
    """
}
