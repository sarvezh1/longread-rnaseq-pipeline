process REFERENCE_VALIDATE {
    tag 'reference compatibility'
    label 'process_single'

    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

    publishDir "${params.outdir}/provenance/reference",
        mode: 'copy',
        pattern: 'reference_validation.*'
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: 'reference_validate.versions.yml'
    publishDir "${params.outdir}/alignments/index",
        mode: params.large_file_publish_mode,
        pattern: '*.fai'

    input:
    path fasta
    path gtf

    output:
    path fasta, emit: validated_fasta
    path "${fasta}.fai", emit: fai
    path 'reference_validation.tsv', emit: report
    path 'reference_validation.sha256', emit: checksums
    path 'reference_validate.versions.yml', emit: versions

    script:
    """
    set -euo pipefail

    samtools faidx ${fasta}

    cut -f1 ${fasta}.fai | sort > fasta.contigs
    if [[ ! -s fasta.contigs ]]; then
        echo 'Reference FASTA contains no indexed sequences.' >&2
        exit 1
    fi
    if [[ \$(uniq -d fasta.contigs | wc -l) -ne 0 ]]; then
        echo 'Reference FASTA contains duplicate sequence names.' >&2
        uniq -d fasta.contigs >&2
        exit 1
    fi

    awk 'BEGIN {FS="\\t"} !/^#/ && NF >= 9 {print \$1}' ${gtf} | sort -u > gtf.contigs
    if [[ ! -s gtf.contigs ]]; then
        echo 'Annotation contains no parseable GTF/GFF records.' >&2
        exit 1
    fi
    comm -23 gtf.contigs fasta.contigs > annotation_only.contigs
    if [[ -s annotation_only.contigs ]]; then
        echo 'Annotation contigs absent from the reference FASTA:' >&2
        sed 's/^/  /' annotation_only.contigs >&2
        exit 1
    fi

    awk 'BEGIN {FS="\\t"; OFS="\\t"}
         NR==FNR {lengths[\$1]=\$2; next}
         !/^#/ && NF >= 9 {
             if (\$4 < 1 || \$5 < \$4 || \$5 > lengths[\$1]) {
                 print \$1, \$4, \$5, lengths[\$1] > "/dev/stderr"
                 bad=1
             }
         }
         END {exit bad}' ${fasta}.fai ${gtf}

    fasta_sha=\$(sha256sum ${fasta} | cut -d' ' -f1)
    gtf_sha=\$(sha256sum ${gtf} | cut -d' ' -f1)
    fai_sha=\$(sha256sum ${fasta}.fai | cut -d' ' -f1)
    fasta_count=\$(wc -l < fasta.contigs)
    gtf_count=\$(wc -l < gtf.contigs)

    printf 'check\\tvalue\\n' > reference_validation.tsv
    printf 'status\\tPASS\\n' >> reference_validation.tsv
    printf 'fasta_file\\t%s\\n' '${fasta}' >> reference_validation.tsv
    printf 'fasta_sha256\\t%s\\n' "\${fasta_sha}" >> reference_validation.tsv
    printf 'fasta_index_sha256\\t%s\\n' "\${fai_sha}" >> reference_validation.tsv
    printf 'annotation_file\\t%s\\n' '${gtf}' >> reference_validation.tsv
    printf 'annotation_sha256\\t%s\\n' "\${gtf_sha}" >> reference_validation.tsv
    printf 'fasta_contigs\\t%s\\n' "\${fasta_count}" >> reference_validation.tsv
    printf 'annotation_contigs\\t%s\\n' "\${gtf_count}" >> reference_validation.tsv
    printf 'annotation_contigs_missing_from_fasta\\t0\\n' >> reference_validation.tsv
    printf '%s  %s\\n%s  %s\\n%s  %s\\n' \
        "\${fasta_sha}" '${fasta}' \
        "\${fai_sha}" '${fasta}.fai' \
        "\${gtf_sha}" '${gtf}' \
        > reference_validation.sha256

    cat > reference_validate.versions.yml <<VERSIONS
    REFERENCE_VALIDATE:
      samtools: "\$(samtools version | sed -n '1s/^samtools //p')"
      sha256sum: "\$(sha256sum --version | sed -n '1s/.* //p')"
      container: "community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd"
    VERSIONS
    """

    stub:
    """
    touch ${fasta}.fai
    printf 'check\\tvalue\\nstatus\\tSTUB\\n' > reference_validation.tsv
    touch reference_validation.sha256
    cat > reference_validate.versions.yml <<VERSIONS
    REFERENCE_VALIDATE:
      samtools: "1.24"
      container: "community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd"
    VERSIONS
    """
}
