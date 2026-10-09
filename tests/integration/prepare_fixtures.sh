#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fixture_dir=${1:-"${root_dir}/.integration-fixtures"}
ont_dir="${fixture_dir}/ont"
pacbio_dir="${fixture_dir}/pacbio"
reference_dir="${fixture_dir}/reference"
mkdir -p "${ont_dir}" "${pacbio_dir}" "${reference_dir}"

samtools_image='community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'
ont_url='https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/071/SRR13619571/SRR13619571_1.fastq.gz'
pacbio_url='https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-SQ2-UHRR-Monomer/1-CCS/m64307e_230628_025302.hifi_reads.bam'
primers_url='https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/REF-primers/IsoSeq_v2_primers_12.fasta'
gencode_url='https://ftp.ebi.ac.uk/pub/databases/gencode/Gencode_human/release_49/gencode.v49.basic.annotation.gtf.gz'
tgfbi_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_000005.10&seq_start=136027988&seq_stop=136064818&rettype=fasta&strand=1'
mitochondrion_url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=nuccore&id=NC_012920.1&rettype=fasta&strand=1'

curl --fail --location --retry 3 --range 0-10485759 \
    --output "${ont_dir}/SRR13619571.prefix10m.fastq.gz" "${ont_url}"
set +o pipefail
gzip -dc "${ont_dir}/SRR13619571.prefix10m.fastq.gz" 2>/dev/null \
    | head -n 16000 > "${ont_dir}/ont.first4000.fastq"
set -o pipefail
[[ $(wc -l < "${ont_dir}/ont.first4000.fastq") -eq 16000 ]]

curl --fail --location --retry 3 --range 0-67108863 \
    --output "${pacbio_dir}/source.prefix64m.bam" "${pacbio_url}"
curl --fail --location --retry 3 \
    --output "${pacbio_dir}/IsoSeq_v2_primers_12.fasta" "${primers_url}"

absolute_fixture_dir=$(cd "${fixture_dir}" && pwd)
docker run --rm -v "${absolute_fixture_dir}:/data" "${samtools_image}" bash -lc '
    set -euo pipefail
    samtools view -H /data/pacbio/source.prefix64m.bam > /tmp/header.sam
    set +o pipefail
    samtools view /data/pacbio/source.prefix64m.bam 2>/tmp/source.stderr \
        | head -n 10000 > /tmp/reads.sam
    set -o pipefail
    cat /tmp/header.sam /tmp/reads.sam \
        | samtools view -b -o /data/pacbio/pacbio.first10000.hifi.bam -
    test "$(samtools view -c /data/pacbio/pacbio.first10000.hifi.bam)" -eq 10000
    samtools quickcheck -u -v /data/pacbio/pacbio.first10000.hifi.bam
'

curl --fail --location --retry 3 --output "${reference_dir}/TGFBI.source.fa" "${tgfbi_url}"
curl --fail --location --retry 3 --output "${reference_dir}/NC_012920.1.source.fa" "${mitochondrion_url}"
curl --fail --location --retry 3 --output "${reference_dir}/gencode.v49.basic.annotation.gtf.gz" "${gencode_url}"

awk 'NR == 1 {print ">chr5_TGFBI"; next} {print}' "${reference_dir}/TGFBI.source.fa" \
    > "${reference_dir}/mini_grch38.fa.lf"
awk 'NR == 1 {print ">chrM"; next} {print}' "${reference_dir}/NC_012920.1.source.fa" \
    >> "${reference_dir}/mini_grch38.fa.lf"
# NCBI records can use different source wrapping and include blank separator
# lines. Rewrap the combined FASTA uniformly for strict readers such as
# pyfaidx (used by IsoQuant), while preserving the sequence exactly.
awk '
    { gsub(/\r/, "") }
    /^>/ {
        if (sequence != "") {
            for (i = 1; i <= length(sequence); i += 60) print substr(sequence, i, 60)
        }
        print
        sequence = ""
        next
    }
    {
        gsub(/[[:space:]]/, "")
        sequence = sequence $0
    }
    END {
        if (sequence != "") {
            for (i = 1; i <= length(sequence); i += 60) print substr(sequence, i, 60)
        }
    }
' "${reference_dir}/mini_grch38.fa.lf" > "${reference_dir}/mini_grch38.fa.rewrapped"
mv "${reference_dir}/mini_grch38.fa.rewrapped" "${reference_dir}/mini_grch38.fa"
rm "${reference_dir}/mini_grch38.fa.lf"

{
    printf '%s\n' '##description: GENCODE v49 TGFBI interval and chrM subset, coordinates localized to mini reference'
    gzip -dc "${reference_dir}/gencode.v49.basic.annotation.gtf.gz" \
        | awk -F '\t' -v OFS='\t' -v start=136027988 '
            !/^#/ && /gene_name "TGFBI"/ {
                $1="chr5_TGFBI"; $4=$4-start+1; $5=$5-start+1; print
            }
            !/^#/ && $1 == "chrM" {print}
        '
} > "${reference_dir}/mini_grch38.gtf.lf"
awk '{printf "%s\r\n", $0}' "${reference_dir}/mini_grch38.gtf.lf" \
    > "${reference_dir}/mini_grch38.gtf"
rm "${reference_dir}/mini_grch38.gtf.lf"

(
    cd "${fixture_dir}"
    sha256sum \
        ont/SRR13619571.prefix10m.fastq.gz \
        ont/ont.first4000.fastq \
        pacbio/source.prefix64m.bam \
        pacbio/pacbio.first10000.hifi.bam \
        pacbio/IsoSeq_v2_primers_12.fasta \
        reference/TGFBI.source.fa \
        reference/NC_012920.1.source.fa \
        reference/gencode.v49.basic.annotation.gtf.gz \
        reference/mini_grch38.fa \
        reference/mini_grch38.gtf \
        > SHA256SUMS
)

cat > "${fixture_dir}/samplesheet.mixed.csv" <<CSV
sample,platform,library_type,preprocessing,condition,replicate,reads,primer_fasta,primer_pair,ont_kit,require_polya
ont_real,ont,pcr_cdna,pychopper,validation,1,${ont_dir}/ont.first4000.fastq,,,PCS109,false
pacbio_real,pacbio,isoseq_hifi,isoseq,validation,1,${pacbio_dir}/pacbio.first10000.hifi.bam,${pacbio_dir}/IsoSeq_v2_primers_12.fasta,IsoSeqX_bc06_5p--IsoSeqX_3p,,true
CSV

cat > "${fixture_dir}/samplesheet.flnc.csv" <<CSV
sample,platform,library_type,preprocessing,condition,replicate,reads,primer_fasta,primer_pair,ont_kit,require_polya
pacbio_flnc_real,pacbio,isoseq_flnc,none,validation,1,${pacbio_dir}/pacbio_real.flnc.bam,,,,false
CSV

printf 'Prepared real integration fixtures in %s\n' "${fixture_dir}"
printf 'Verify recorded checksums with: (cd %s && sha256sum -c SHA256SUMS)\n' "${fixture_dir}"
