process NORMALIZED_MANIFEST {
    tag "${meta.id} (normalized aligned reads)"
    label 'process_single'

    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

    publishDir "${params.outdir}/provenance/samples",
        mode: 'copy',
        pattern: "*.normalized.{json,tsv}"
    publishDir "${params.outdir}/provenance/samples",
        mode: 'copy',
        pattern: "*.runtime_versions.yml"

    input:
    tuple val(meta), path(bam), path(bai), path(qc_dir), path(version_files, stageAs: 'versions/*')
    path reference_validation
    val fasta_source
    val gtf_source
    val genome_build
    val annotation_release

    output:
    tuple val(meta), path(bam), path(bai), path(qc_dir), path("${meta.id}.normalized.json"), emit: normalized
    tuple val(meta), path("${meta.id}.normalized.tsv"), emit: manifest_tsv
    tuple val(meta), path("${meta.id}.runtime_versions.yml"), emit: runtime_versions

    script:
    def preprocessing_steps = meta.platform == 'ont' ?
        (meta.preprocessing == 'pychopper' ? ['NanoPlot raw QC', 'Pychopper', 'minimap2 splice alignment'] : ['NanoPlot raw QC', 'Pychopper bypassed', 'minimap2 splice alignment']) :
        (meta.library_type == 'isoseq_hifi' ? ['NanoPlot raw QC', 'lima --isoseq', 'isoseq refine', 'isoseq cluster2', 'pbmm2 ISOSEQ alignment'] : ['NanoPlot raw QC', 'lima/refine bypassed (FLNC input)', 'isoseq cluster2', 'pbmm2 ISOSEQ alignment'])
    def alignment_command = meta.platform == 'ont' ?
        "minimap2 -ax splice -R <read-group> <splice-index> <reads> | samtools sort; samtools index" :
        "pbmm2 align <ISOSEQ-index> <clustered-bam> <output-bam> --preset ISOSEQ --unmapped --sort --bam-index BAI"
    def manifest = [
        contract_version: '1.0',
        sample: meta.id,
        platform: meta.platform,
        library_type: meta.library_type,
        condition: meta.condition,
        replicate: meta.replicate,
        original_input: meta.original_input,
        sorted_bam: "alignments/${meta.platform}/${bam.name}",
        bam_index: "alignments/${meta.platform}/${bai.name}",
        alignment_qc_directory: "qc/alignment/${qc_dir.name}",
        reference: [
            genome_build: genome_build,
            fasta: fasta_source,
            annotation_gtf: gtf_source,
            annotation_release: annotation_release
        ],
        preprocessing: [
            selection: meta.preprocessing,
            steps: preprocessing_steps,
            ont_kit: meta.ont_kit ?: null,
            require_polya: meta.require_polya,
            include_singletons: params.isoseq_include_singletons,
            primer_pair: meta.primer_pair ?: null
        ],
        alignment: [
            command_template: alignment_command,
            annotation_guided_alignment: false
        ],
        runtime_tool_versions: "provenance/samples/${meta.id}.runtime_versions.yml",
        reference_validation: "provenance/reference/${reference_validation.name}"
    ]
    def manifest_json = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson(manifest))
    def manifest_tsv = [
        'sample\tplatform\tlibrary_type\tcondition\treplicate\tsorted_bam\tbam_index\talignment_qc\tgenome_build\tannotation_release',
        [
            meta.id,
            meta.platform,
            meta.library_type,
            meta.condition,
            meta.replicate,
            "alignments/${meta.platform}/${bam.name}",
            "alignments/${meta.platform}/${bai.name}",
            "qc/alignment/${qc_dir.name}",
            genome_build,
            annotation_release
        ].join('\t')
    ].join('\n')
    def manifest_json_b64 = manifest_json.bytes.encodeBase64().toString()
    def manifest_tsv_b64 = (manifest_tsv + '\n').bytes.encodeBase64().toString()
    """
    cat versions/*.yml > ${meta.id}.runtime_versions.yml
    printf '%s' '${manifest_json_b64}' | base64 --decode > ${meta.id}.normalized.json
    printf '%s' '${manifest_tsv_b64}' | base64 --decode > ${meta.id}.normalized.tsv
    """

    stub:
    def stub_manifest = groovy.json.JsonOutput.prettyPrint(groovy.json.JsonOutput.toJson([
        contract_version: '1.0',
        sample: meta.id,
        platform: meta.platform,
        library_type: meta.library_type,
        sorted_bam: bam.name,
        bam_index: bai.name
    ]))
    def stub_manifest_json_b64 = stub_manifest.bytes.encodeBase64().toString()
    def stub_manifest_tsv = "sample\tplatform\tlibrary_type\n${meta.id}\t${meta.platform}\t${meta.library_type}\n"
    def stub_manifest_tsv_b64 = stub_manifest_tsv.bytes.encodeBase64().toString()
    """
    touch ${meta.id}.runtime_versions.yml
    printf '%s' '${stub_manifest_json_b64}' | base64 --decode > ${meta.id}.normalized.json
    printf '%s' '${stub_manifest_tsv_b64}' | base64 --decode > ${meta.id}.normalized.tsv
    """
}
