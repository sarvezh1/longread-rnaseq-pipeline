#!/usr/bin/env bash
set -euo pipefail

script_dir=$(CDPATH= cd -- "$(dirname -- "$0")" && pwd)
project_root=$(CDPATH= cd -- "${script_dir}/../.." && pwd)
test_root=$(mktemp -d)
trap 'rm -rf "${test_root}"' EXIT
cd "${project_root}"

for caller in isoquant flair bambu; do
    python3 bin/standardize_benchmark_transcripts.py \
        --sample fixture --platform ont --condition unit --replicate 1 --caller "${caller}" \
        --gtf "tests/data/benchmark/${caller}.gtf" \
        --classification "tests/data/benchmark/${caller}.classification.tsv" \
        --support tests/data/benchmark/support.tsv \
        --support-kind count \
        --transcripts "${test_root}/${caller}.tsv" \
        --junctions "${test_root}/${caller}.junctions.tsv" \
        --summary "${test_root}/${caller}.json"
done

python3 bin/compare_transcript_structures.py \
    --sample fixture --platform ont --condition unit --replicate 1 \
    --isoquant "${test_root}/isoquant.tsv" --flair "${test_root}/flair.tsv" --bambu "${test_root}/bambu.tsv" \
    --membership "${test_root}/membership.tsv" --pairwise "${test_root}/pairwise.tsv" \
    --ends "${test_root}/ends.tsv" --partial-junctions "${test_root}/partial.tsv" \
    --gene-overlap "${test_root}/genes.tsv" --summary "${test_root}/summary.json"

python3 - "${test_root}" <<'PY'
import csv, json, pathlib, sys
root = pathlib.Path(sys.argv[1])
summary = json.loads((root / "summary.json").read_text())
assert summary["membership_counts"]["isoquant+flair+bambu"] == 2  # exact multi-exon + single-exon overlap
assert summary["membership_counts"]["isoquant_only"] >= 1
assert summary["membership_counts"]["flair_only"] >= 1
assert summary["membership_counts"]["bambu_only"] >= 1
assert summary["pairwise"]["isoquant_vs_flair"]["exact_multi_exon_structures"]["intersection"] == 1
assert summary["pairwise"]["isoquant_vs_flair"]["intron_chains"]["intersection"] == 2
ends = list(csv.DictReader((root / "ends.tsv").open(), delimiter="\t"))
chain = next(row for row in ends if row["left_transcript_id"] == "iq_chain" and row["right_transcript_id"] == "fl_chain")
assert chain["tss_difference"] == "-10" and chain["tes_difference"] == "21"
partial = list(csv.DictReader((root / "partial.tsv").open(), delimiter="\t"))
assert any(row["left_transcript_id"] == "iq_partial" and row["right_transcript_id"] == "fl_partial" for row in partial)
assert summary["consensus_transcriptome_created"] is False
PY

# Empty outputs must be valid and all empty-set Jaccards must be explicit NA/null.
head -n 1 "${test_root}/isoquant.tsv" > "${test_root}/empty.tsv"
python3 bin/compare_transcript_structures.py \
    --sample empty --platform pacbio --condition unit --replicate 1 \
    --isoquant "${test_root}/empty.tsv" --flair "${test_root}/empty.tsv" --bambu "${test_root}/empty.tsv" \
    --membership "${test_root}/empty.membership.tsv" --pairwise "${test_root}/empty.pairwise.tsv" \
    --ends "${test_root}/empty.ends.tsv" --partial-junctions "${test_root}/empty.partial.tsv" \
    --gene-overlap "${test_root}/empty.genes.tsv" --summary "${test_root}/empty.json"
python3 - "${test_root}/empty.json" <<'PY'
import json, sys
data = json.load(open(sys.argv[1]))
for pair in data["pairwise"].values():
    for metric in pair.values():
        assert metric["union"] == 0 and metric["jaccard"] is None and metric["status"] == "both_empty"
PY

# DSL2 routing: one normalized record must fan out to the reused IsoQuant,
# FLAIR, and Bambu modules, then through three common SQANTI3 instances.
nextflow run . -profile test_stub -stub-run \
    -w "${test_root}/benchmark.work" \
    --outdir "${test_root}/benchmark.results" \
    --mode benchmark \
    --input tests/data/structural/ont_valid.csv
for caller in isoquant flair bambu; do
    test -d "${test_root}/benchmark.results/benchmark/${caller}"
    test -d "${test_root}/benchmark.results/benchmark/sqanti3/${caller}"
    test -s "${test_root}/benchmark.results/benchmark/concordance/transcripts/ont_only_stub.${caller}.standardized_transcripts.tsv"
done
test -s "${test_root}/benchmark.results/benchmark/summary/ont_only_stub.benchmark_contract.json"
test -s "${test_root}/benchmark.results/provenance/software_versions/ont_only_stub.benchmark.versions.yml"

echo 'Block 3 structural concordance contracts passed.'
