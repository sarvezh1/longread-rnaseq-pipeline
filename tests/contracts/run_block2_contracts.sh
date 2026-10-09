#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "${script_dir}/../.." && pwd)
test_root=$(mktemp -d)
trap 'rm -rf "${test_root}"' EXIT

cd "${project_root}"

python3 bin/summarize_transcript_structure.py \
    --gtf tests/data/structural/block2_valid.gtf \
    --classification tests/data/structural/block2_classification.tsv \
    --support-counts tests/data/structural/block2_counts.tsv \
    --sample structural_unit \
    --structure-output "${test_root}/structure.tsv" \
    --json-output "${test_root}/summary.json" \
    --markdown-output "${test_root}/summary.md"

grep -q $'tx_single\tgene2\tchr1\t-\t701\t800\t100\t1\tsingle-exon' "${test_root}/structure.tsv"
grep -q '"novel_in_catalog": 1' "${test_root}/summary.json"
grep -q '"single_exon_transcripts": 1' "${test_root}/summary.json"

if python3 bin/summarize_transcript_structure.py \
    --gtf tests/data/structural/block2_malformed.gtf \
    --classification tests/data/structural/block2_classification.tsv \
    --support-counts tests/data/structural/block2_counts.tsv \
    --sample malformed \
    --structure-output "${test_root}/bad.tsv" \
    --json-output "${test_root}/bad.json" \
    --markdown-output "${test_root}/bad.md" 2>/dev/null; then
    echo 'Malformed transcript GTF unexpectedly passed validation.' >&2
    exit 1
fi

cp tests/data/structural/reference.fa "${test_root}/reference.fa"
cp tests/data/structural/block2_valid.gtf "${test_root}/reference.gtf"
fasta_sha=$(sha256sum "${test_root}/reference.fa" | cut -d' ' -f1)
gtf_sha=$(sha256sum "${test_root}/reference.gtf" | cut -d' ' -f1)
cat > "${test_root}/reference_validation.tsv" <<EOF
check	value
status	PASS
fasta_sha256	${fasta_sha}
annotation_sha256	${gtf_sha}
EOF
cat > "${test_root}/normalized.json" <<'EOF'
{
  "contract_version": "1.0",
  "sample": "contract_unit",
  "platform": "ont",
  "library_type": "cdna",
  "condition": "unit",
  "replicate": "1",
  "original_input": "reads.fastq",
  "sorted_bam": "alignments/ont/sample.bam",
  "bam_index": "alignments/ont/sample.bam.bai",
  "alignment_qc_directory": "qc/alignment/sample",
  "reference": {"genome_build": "unit-build", "annotation_release": "unit-release"},
  "preprocessing": {},
  "alignment": {},
  "runtime_tool_versions": "versions.yml",
  "reference_validation": "reference_validation.tsv"
}
EOF

python3 bin/validate_block2_contract.py \
    --manifest "${test_root}/normalized.json" \
    --reference-validation "${test_root}/reference_validation.tsv" \
    --fasta "${test_root}/reference.fa" \
    --gtf "${test_root}/reference.gtf" \
    --sample contract_unit \
    --platform ont \
    --library-type cdna \
    --genome-build unit-build \
    --annotation-release unit-release \
    --output "${test_root}/validation.json"
grep -q '"status": "PASS"' "${test_root}/validation.json"

if python3 bin/validate_block2_contract.py \
    --manifest "${test_root}/normalized.json" \
    --reference-validation "${test_root}/reference_validation.tsv" \
    --fasta "${test_root}/reference.fa" \
    --gtf "${test_root}/reference.gtf" \
    --sample contract_unit \
    --platform ont \
    --library-type cdna \
    --genome-build wrong-build \
    --annotation-release unit-release \
    --output "${test_root}/mismatch.json" 2>/dev/null; then
    echo 'Incompatible reference metadata unexpectedly passed validation.' >&2
    exit 1
fi

python3 - "${test_root}/normalized.json" "${test_root}/missing.json" <<'PY'
import json
import sys
data = json.load(open(sys.argv[1]))
del data['bam_index']
json.dump(data, open(sys.argv[2], 'w'))
PY
if python3 bin/validate_block2_contract.py \
    --manifest "${test_root}/missing.json" \
    --reference-validation "${test_root}/reference_validation.tsv" \
    --fasta "${test_root}/reference.fa" \
    --gtf "${test_root}/reference.gtf" \
    --sample contract_unit \
    --platform ont \
    --library-type cdna \
    --genome-build unit-build \
    --annotation-release unit-release \
    --output "${test_root}/missing-validation.json" 2>/dev/null; then
    echo 'Incomplete Block 1 contract unexpectedly passed validation.' >&2
    exit 1
fi

python3 bin/build_transcriptome_contract.py \
    --sample contract_unit \
    --platform ont \
    --library-type cdna \
    --condition unit \
    --replicate 1 \
    --genome-build unit-build \
    --annotation-release unit-release \
    --fasta-sha256 "${fasta_sha}" \
    --annotation-sha256 "${gtf_sha}" \
    --isoquant-data-type nanopore \
    --isoquant-full-length false \
    --isoquant-novel-unspliced tool-native \
    --primary-gtf transcriptome/raw/contract_unit.raw.gtf \
    --transcript-fasta transcriptome/raw/contract_unit.raw.fasta \
    --classification isoforms/sqanti3/classification.tsv \
    --junctions isoforms/sqanti3/junctions.tsv \
    --support-table isoforms/isoquant/support.tsv.gz \
    --counts-directory isoforms/isoquant \
    --sqanti-directory isoforms/sqanti3 \
    --structure-table qc/transcript_structure/structure.tsv \
    --summary-json qc/transcript_structure/summary.json \
    --block1-manifest provenance/samples/normalized.json \
    --isoquant-versions provenance/software_versions/isoquant.yml \
    --sqanti3-versions provenance/software_versions/sqanti3.yml \
    --output "${test_root}/transcriptome_contract.json"
grep -q '"isoquant_versions": "provenance/software_versions/isoquant.yml"' "${test_root}/transcriptome_contract.json"
grep -q '"sqanti3_versions": "provenance/software_versions/sqanti3.yml"' "${test_root}/transcriptome_contract.json"
grep -q '"curated_transcriptome": null' "${test_root}/transcriptome_contract.json"

echo 'Block 2 contract and transcript-structure tests passed.'
