process SAMTOOLS_ALIGNMENT_QC {
    tag "${meta.id} (${meta.platform})"
    label 'process_medium'

    container 'community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

    publishDir "${params.outdir}/qc/alignment",
        mode: 'copy',
        pattern: "*.alignment_qc"
    publishDir "${params.outdir}/provenance/software_versions",
        mode: 'copy',
        pattern: "*.samtools_alignment_qc.versions.yml"

    input:
    tuple val(meta), path(bam), path(bai)
    path fasta
    path fai

    output:
    tuple val(meta), path("${meta.id}.alignment_qc"), emit: qc
    tuple val(meta), path("${meta.id}.samtools_alignment_qc.versions.yml"), emit: versions

    script:
    """
    mkdir -p ${meta.id}.alignment_qc

    test -s ${bam}
    test -s ${bai}
    samtools quickcheck -v ${bam}
    sort_order=\$(samtools view -H ${bam} | awk -F'\\t' '/^@HD/ {for (i=1; i<=NF; i++) if (\$i ~ /^SO:/) {sub(/^SO:/, "", \$i); print \$i; exit}}')
    if [[ "\${sort_order}" != 'coordinate' ]]; then
        echo "BAM for sample ${meta.id} is not declared coordinate sorted (SO=\${sort_order:-missing})." >&2
        exit 1
    fi

    samtools view -H ${bam} | awk -F'\\t' '/^@SQ/ {for (i=1; i<=NF; i++) if (\$i ~ /^SN:/) {sub(/^SN:/, "", \$i); print \$i}}' | sort -u > bam.contigs
    cut -f1 ${fai} | sort -u > reference.contigs
    comm -23 bam.contigs reference.contigs > bam_only.contigs
    comm -13 bam.contigs reference.contigs > reference_only.contigs
    if [[ -s bam_only.contigs || -s reference_only.contigs ]]; then
        echo "BAM/reference sequence dictionaries differ for sample ${meta.id}." >&2
        sed 's/^/BAM only: /' bam_only.contigs >&2
        sed 's/^/Reference only: /' reference_only.contigs >&2
        exit 1
    fi

    expected_platform='${meta.platform == 'ont' ? 'ONT' : 'PACBIO'}'
    if ! samtools view -H ${bam} | awk -F'\\t' -v sample='${meta.id}' -v platform="\${expected_platform}" '
        /^@RG/ {
            sm=""; pl="";
            for (i=1; i<=NF; i++) {
                if (\$i ~ /^SM:/) {sm=\$i; sub(/^SM:/, "", sm)}
                if (\$i ~ /^PL:/) {pl=\$i; sub(/^PL:/, "", pl)}
            }
            if (sm == sample && pl == platform) found=1
        }
        END {exit(found ? 0 : 1)}'; then
        echo "BAM for sample ${meta.id} lacks the expected SM=${meta.id}, PL=\${expected_platform} read group." >&2
        exit 1
    fi

    samtools flagstat -@ ${task.cpus} ${bam} > ${meta.id}.alignment_qc/${meta.id}.flagstat.txt
    samtools stats --threads ${task.cpus} --reference ${fasta} ${bam} > ${meta.id}.alignment_qc/${meta.id}.stats.txt
    samtools idxstats ${bam} > ${meta.id}.alignment_qc/${meta.id}.idxstats.tsv
    samtools coverage ${bam} > ${meta.id}.alignment_qc/${meta.id}.coverage.tsv

    total=\$(samtools view -c ${bam})
    if [[ "\${total}" -eq 0 ]]; then
        echo "BAM for sample ${meta.id} contains no records." >&2
        exit 1
    fi
    mapped=\$(samtools view -c -F 4 ${bam})
    primary=\$(samtools view -c -F 2308 ${bam})
    secondary=\$(samtools view -c -f 256 ${bam})
    supplementary=\$(samtools view -c -f 2048 ${bam})
    unmapped=\$(samtools view -c -f 4 ${bam})
    mapped_pct=\$(awk -v m="\${mapped}" -v t="\${total}" 'BEGIN { if (t == 0) print "0.0000"; else printf "%.4f", (100*m/t) }')
    mean_mapq=\$(samtools view -F 4 ${bam} | awk '{sum += \$5; n += 1} END {if (n == 0) print "0.0000"; else printf "%.4f", sum/n}')
    primary_total=\$(samtools view -c -F 2304 ${bam})
    rg_count=\$(samtools view -H ${bam} | grep -c '^@RG')

    printf 'metric\\tvalue\\n' > ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'total_records\\t%s\\n' "\${total}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'mapped_records\\t%s\\n' "\${mapped}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'mapped_percentage\\t%s\\n' "\${mapped_pct}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'primary_mapped_alignments\\t%s\\n' "\${primary}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'primary_records\\t%s\\n' "\${primary_total}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'secondary_alignments\\t%s\\n' "\${secondary}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'supplementary_alignments\\t%s\\n' "\${supplementary}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'unmapped_records\\t%s\\n' "\${unmapped}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    printf 'mean_mapping_quality_mapped_records\\t%s\\n' "\${mean_mapq}" >> ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv

    printf 'check\\tvalue\\n' > ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'status\\tPASS\\n' >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'quickcheck\\tPASS\\n' >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'coordinate_sort_order\\t%s\\n' "\${sort_order}" >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'bam_index_readable\\tPASS\\n' >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'reference_dictionary_match\\tPASS\\n' >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'expected_sample_read_group\\tPASS\\n' >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'read_group_count\\t%s\\n' "\${rg_count}" >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    printf 'non_empty_records\\t%s\\n' "\${total}" >> ${meta.id}.alignment_qc/${meta.id}.validation.tsv

    printf 'mapping_quality\\trecords\\n' > ${meta.id}.alignment_qc/${meta.id}.mapq_histogram.tsv
    samtools view -F 4 ${bam} \
        | awk '{counts[\$5] += 1} END {for (q in counts) print q "\t" counts[q]}' \
        | sort -n -k1,1 \
        >> ${meta.id}.alignment_qc/${meta.id}.mapq_histogram.tsv

    cat > ${meta.id}.samtools_alignment_qc.versions.yml <<VERSIONS
    SAMTOOLS_ALIGNMENT_QC:
      samtools: "\$(samtools version | sed -n '1s/^samtools //p')"
      acceptance_thresholds: null
      container: "community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd"
    VERSIONS
    """

    stub:
    """
    mkdir -p ${meta.id}.alignment_qc
    printf 'metric\\tvalue\\n' > ${meta.id}.alignment_qc/${meta.id}.alignment_summary.tsv
    touch ${meta.id}.alignment_qc/${meta.id}.flagstat.txt
    touch ${meta.id}.alignment_qc/${meta.id}.stats.txt
    touch ${meta.id}.alignment_qc/${meta.id}.idxstats.tsv
    touch ${meta.id}.alignment_qc/${meta.id}.coverage.tsv
    touch ${meta.id}.alignment_qc/${meta.id}.mapq_histogram.tsv
    printf 'check\\tvalue\\nstatus\\tSTUB\\n' > ${meta.id}.alignment_qc/${meta.id}.validation.tsv
    cat > ${meta.id}.samtools_alignment_qc.versions.yml <<VERSIONS
    SAMTOOLS_ALIGNMENT_QC:
      samtools: "1.24"
      acceptance_thresholds: null
      container: "community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd"
    VERSIONS
    """
}
