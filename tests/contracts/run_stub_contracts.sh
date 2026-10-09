#!/usr/bin/env bash

set -euo pipefail

script_dir="$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)"
project_root="$(CDPATH= cd -- "${script_dir}/../.." && pwd)"
test_root="$(mktemp -d)"
trap 'rm -rf "${test_root}"' EXIT

cd "${project_root}"

python3 tests/contracts/audit_fixture_paths.py

if ! command -v nextflow >/dev/null 2>&1; then
        echo "Nextflow is required to run the Block 1 stub contracts." >&2
    exit 127
fi

run_expected_failure() {
    case_name="$1"
    shift
    if nextflow run . -profile test_stub -stub-run -w "${test_root}/${case_name}.work" --outdir "${test_root}/${case_name}.results" "$@"; then
        echo "Expected failure did not occur: ${case_name}" >&2
        exit 1
    fi
}

nextflow run . \
    -profile test_stub \
    -stub-run \
    -w "${test_root}/ont-only.work" \
    --outdir "${test_root}/ont-only.results" \
    --input tests/data/structural/ont_valid.csv
test -s "${test_root}/ont-only.results/provenance/samples/ont_only_stub.transcriptome_contract.json"
test -s "${test_root}/ont-only.results/qc/transcript_structure/ont_only_stub.transcript_summary.json"

nextflow run . \
    -profile test_stub \
    -stub-run \
    -w "${test_root}/pacbio-only.work" \
    --outdir "${test_root}/pacbio-only.results" \
    --input tests/data/structural/pacbio_valid.csv
test -s "${test_root}/pacbio-only.results/provenance/samples/pacbio_only_stub.transcriptome_contract.json"

nextflow run . \
    -profile test_stub \
    -stub-run \
    -w "${test_root}/mixed.work" \
    --outdir "${test_root}/mixed.results"

nextflow run . \
    -profile test_stub \
    -stub-run \
    -w "${test_root}/bypass.work" \
    --outdir "${test_root}/bypass.results" \
    --input tests/data/structural/bypass_valid.csv

run_expected_failure unsupported_combination \
    --input tests/data/structural/unsupported_combination.csv

run_expected_failure missing_primer \
    --input tests/data/structural/missing_primer.csv

run_expected_failure missing_reference \
    --fasta "${test_root}/does-not-exist.fa"

echo "Block 1 and Block 2 stub contracts passed."
