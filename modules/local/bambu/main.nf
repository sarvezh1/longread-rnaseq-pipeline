process BAMBU {
    tag "${meta.id} (${meta.platform})"
    label 'process_high'

    container 'longread-rnaseq-pipeline/bambu:3.14.0-bioc3.23-r4.6'

    publishDir "${params.outdir}/benchmark/bambu", mode: 'copy', pattern: '*.bambu'
    publishDir "${params.outdir}/provenance/software_versions", mode: 'copy', pattern: '*.bambu.versions.yml'

    input:
    tuple val(meta), path(bam), path(bai), path(alignment_qc), path(block1_manifest), path(input_validation)
    path fasta
    path gtf

    output:
    tuple val(meta), path("${meta.id}.bambu"), path("${meta.id}.bambu.gtf"), path("${meta.id}.bambu.support.tsv"), path("${meta.id}.bambu.counts.tsv"), path(block1_manifest), path(input_validation), emit: results
    tuple val(meta), path("${meta.id}.bambu.versions.yml"), emit: versions

    script:
    """
    set -euo pipefail
    export NXF_TASK_CPUS=${task.cpus}
    mkdir -p ${meta.id}.bambu
    Rscript ${projectDir}/bin/run_bambu.R ${bam} ${bai} ${gtf} ${fasta} ${meta.id}.bambu 2>&1 | tee ${meta.id}.bambu/execution.log
    detected_gtf='${meta.id}.bambu/detected_models/detected_extended_annotations.gtf'
    transcript_counts='${meta.id}.bambu/detected_models/bambu_transcript_counts.tsv'
    test -n "\${detected_gtf}" && test -f "\${detected_gtf}"
    test -n "\${transcript_counts}" && test -f "\${transcript_counts}"
    cp "\${detected_gtf}" ${meta.id}.bambu.gtf
    cp "\${transcript_counts}" ${meta.id}.bambu.counts.tsv
    cp "\${transcript_counts}" ${meta.id}.bambu.support.tsv
    bambu_version=\$(Rscript -e 'cat(as.character(packageVersion("bambu")))')
    r_version=\$(Rscript -e 'cat(paste(R.version\$major, R.version\$minor, sep="."))')
    bioc_version=\$(Rscript -e 'cat(as.character(BiocManager::version()))')
    quantification_status=\$(awk -F '\t' '\$1 == "status" {print \$2; exit}' ${meta.id}.bambu/quantification_status.tsv)
    printf '%s\n' \
        'BAMBU:' \
        "  bambu: \"\${bambu_version}\"" \
        "  bioconductor: \"\${bioc_version}\"" \
        "  r: \"\${r_version}\"" \
        '  container: "longread-rnaseq-pipeline/bambu:3.14.0-bioc3.23-r4.6"' \
        '  base_container: "bioconductor/bioconductor_docker:RELEASE_3_23@sha256:821dbf9ac119eac41f177531c7ca8fc7084c99c04eb218a1aa7f52cb63bad5d9"' \
        '  discovery: true' \
        '  quantification: true' \
        "  quantification_status: \"\${quantification_status}\"" \
        '  comparison_universe: "Bambu-native readCount >= 1"' \
        > ${meta.id}.bambu.versions.yml
    """

    stub:
    """
    mkdir -p ${meta.id}.bambu/detected_models ${meta.id}.bambu/all_models
    printf 'chrStub\tBambu\ttranscript\t1\t100\t.\t+\t.\tgene_id "bambu_gene"; transcript_id "bambu_tx";\nchrStub\tBambu\texon\t1\t100\t.\t+\t.\tgene_id "bambu_gene"; transcript_id "bambu_tx";\n' > ${meta.id}.bambu.gtf
    printf 'transcript_id\tcount\nbambu_tx\t1\n' > ${meta.id}.bambu.counts.tsv
    cp ${meta.id}.bambu.counts.tsv ${meta.id}.bambu.support.tsv
    cp ${meta.id}.bambu.gtf ${meta.id}.bambu/detected_models/
    cp ${meta.id}.bambu.gtf ${meta.id}.bambu/all_models/
    printf 'BAMBU:\n  bambu: "3.14.0"\n  bioconductor: "3.23"\n  r: "4.6"\n' > ${meta.id}.bambu.versions.yml
    """
}
