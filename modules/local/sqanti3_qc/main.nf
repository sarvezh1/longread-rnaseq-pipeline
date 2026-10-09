process SQANTI3_QC {
    tag "${meta.id} (structural QC)"
    label 'process_high'

    container 'quay.io/biocontainers/sqanti3:6.0.2--hdfd78af_1@sha256:3ebeb658fb2e0e308e73cba50dd3b97c63963d8225f303f298ad8ee76264939d'

    publishDir { "${params.outdir}/${meta.benchmark_caller ? 'benchmark/sqanti3/' + meta.benchmark_caller : 'isoforms/sqanti3/per_sample'}" },
        mode: 'copy',
        pattern: '*.sqanti3'
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: '*.sqanti3.versions.yml'

    input:
    tuple val(meta), path(isoquant_dir), path(primary_gtf), path(support), path(counts), path(block1_manifest), path(input_validation)
    path fasta
    path gtf

    output:
    tuple val(meta), path(isoquant_dir), path(primary_gtf), path(support), path(counts), path("${meta.id}.sqanti3"), path("${meta.id}.sqanti3_classification.txt"), path("${meta.id}.sqanti3_junctions.txt"), path("${meta.id}.transcripts.fasta"), path(block1_manifest), path(input_validation), emit: results
    tuple val(meta), path("${meta.id}.sqanti3.versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    mkdir -p ${meta.id}.sqanti3
    mkdir -p sqanti3.compat
    ln -s /usr/local/lib/libbz2.so.1.0 sqanti3.compat/libbz2.so.1
    export LD_LIBRARY_PATH="\$PWD/sqanti3.compat:/usr/local/lib:\${LD_LIBRARY_PATH:-}"

    sqanti_version=\$(sqanti3_qc.py -v 2>&1 | awk 'NF { version=\$NF } END { print version }')
    if awk 'BEGIN { FS="\\t"; found=0 } !/^#/ && (\$3 == "transcript" || \$3 == "exon") { found=1; exit } END { exit !found }' ${primary_gtf}; then
        sqanti3_qc.py \
            --isoforms ${primary_gtf} \
            --refGTF ${gtf} \
            --refFasta ${fasta} \
            --report skip \
            -t ${task.cpus} \
            -n 1 \
            -o ${meta.id}.sqanti3 \
            -d ${meta.id}.sqanti3

        test -s ${meta.id}.sqanti3/${meta.id}.sqanti3_classification.txt
        test -s ${meta.id}.sqanti3/${meta.id}.sqanti3_junctions.txt
        test -s ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.gtf
        test -s ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.fasta
        cp ${meta.id}.sqanti3/${meta.id}.sqanti3_classification.txt ${meta.id}.sqanti3_classification.txt
        cp ${meta.id}.sqanti3/${meta.id}.sqanti3_junctions.txt ${meta.id}.sqanti3_junctions.txt
        cp ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.fasta ${meta.id}.transcripts.fasta
        sqanti_status='CLASSIFIED'
    else
        printf 'isoform\tchrom\tstrand\tlength\texons\tstructural_category\tassociated_gene\tassociated_transcript\n' > ${meta.id}.sqanti3_classification.txt
        printf 'isoform\tchrom\tstrand\n' > ${meta.id}.sqanti3_junctions.txt
        : > ${meta.id}.transcripts.fasta
        cp ${primary_gtf} ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.gtf
        cp ${meta.id}.transcripts.fasta ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.fasta
        cp ${meta.id}.sqanti3_classification.txt ${meta.id}.sqanti3/
        cp ${meta.id}.sqanti3_junctions.txt ${meta.id}.sqanti3/
        printf 'NO_TRANSCRIPT_MODELS\n' > ${meta.id}.sqanti3/${meta.id}.sqanti3_status.txt
        sqanti_status='NO_TRANSCRIPT_MODELS'
    fi

    {
    printf '%s\\n' \
        'SQANTI3_QC:'
    printf '  sqanti3: "%s"\\n' "\${sqanti_version}"
    printf '%s\\n' \
        '  container: "quay.io/biocontainers/sqanti3:6.0.2--hdfd78af_1@sha256:3ebeb658fb2e0e308e73cba50dd3b97c63963d8225f303f298ad8ee76264939d"' \
        '  classification_schema: "6.x"' \
        '  input_type: "GTF"' \
        '  orf_prediction: "skipped"' \
        '  report: "skipped"' \
        '  container_compatibility: "task-local libbz2.so.1 symlink to /usr/local/lib/libbz2.so.1.0"'
    printf '  status: "%s"\\n' "\${sqanti_status}"
    } > ${meta.id}.sqanti3.versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}.sqanti3
    printf '%s\n' \
        'isoform\tchrom\tstrand\tlength\texons\tstructural_category\tassociated_gene\tassociated_transcript' \
        'stub_transcript\tchrStub\t+\t100\t1\tfull-splice_match\tstub_gene\tstub_reference' \
        > ${meta.id}.sqanti3_classification.txt
    printf 'isoform\tchrom\tstrand\n' > ${meta.id}.sqanti3_junctions.txt
    printf '>stub_transcript\nACGT\n' > ${meta.id}.transcripts.fasta
    cp ${meta.id}.sqanti3_classification.txt ${meta.id}.sqanti3/
    cp ${meta.id}.sqanti3_junctions.txt ${meta.id}.sqanti3/
    cp ${meta.id}.transcripts.fasta ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.fasta
    cp ${primary_gtf} ${meta.id}.sqanti3/${meta.id}.sqanti3_corrected.gtf
    printf '%s\n' \
        'SQANTI3_QC:' \
        '  sqanti3: "6.0.2"' \
        '  container: "quay.io/biocontainers/sqanti3:6.0.2--hdfd78af_1@sha256:3ebeb658fb2e0e308e73cba50dd3b97c63963d8225f303f298ad8ee76264939d"' \
        '  classification_schema: "6.x"' \
        '  orf_prediction: "skipped"' \
        '  status: "STUB"' \
        > ${meta.id}.sqanti3.versions.yml
    """
}
