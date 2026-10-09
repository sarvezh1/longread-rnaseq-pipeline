#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
fixture_dir=${1:-"${root_dir}/.integration-fixtures"}
out_dir=${2:-"${fixture_dir}/results"}
work_dir=${3:-"${fixture_dir}/work"}
pipeline_mode=${4:-default}

test -s "${fixture_dir}/samplesheet.mixed.csv"
test -s "${fixture_dir}/reference/mini_grch38.fa"
test -s "${fixture_dir}/reference/mini_grch38.gtf"

cd "${root_dir}"
nextflow run . \
    -profile docker \
    -c tests/integration/nextflow.config \
    -w "${work_dir}" \
    --input "${fixture_dir}/samplesheet.mixed.csv" \
    --outdir "${out_dir}" \
    --fasta "${fixture_dir}/reference/mini_grch38.fa" \
    --gtf "${fixture_dir}/reference/mini_grch38.gtf" \
    --genome_build 'GRCh38.p14-mini-validation' \
    --annotation_release 'GENCODE 49 subset' \
    --lima_peek_guess true \
    --mode "${pipeline_mode}"
