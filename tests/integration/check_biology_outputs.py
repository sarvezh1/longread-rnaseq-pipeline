#!/usr/bin/env python3
"""Validate the real compact-fixture Block 4 contract without biological inference."""
import csv, json, sys
from pathlib import Path

root=Path(sys.argv[1] if len(sys.argv)>1 else ".integration-fixtures/block4-results")
contract=json.loads((root/"biology/biology_contract.json").read_text())
assert contract["production_transcriptome"].startswith("IsoQuant")
assert contract["benchmark_callers_used_to_define_transcriptome"] is False
for name,relative in contract["paths"].items():
    assert (root/relative).exists(), f"missing contract path {name}: {relative}"

status=contract["analysis_status"]
assert status["design"]["status"]=="NOT_RUN_NO_EXPLICIT_CONTRAST"
for key in ("gene_differential_expression","transcript_differential_expression","differential_transcript_usage","isoform_switches"):
    assert status[key]["status"]=="NOT_RUN_NO_EXPLICIT_CONTRAST"
assert status["differential_splicing"]["status"].startswith("NOT_RUN_")
assert status["expression_qc"]["sample_exclusion_performed"] is False
assert status["orfs"]["prediction_status"]=="PREDICTION_NOT_AVAILABLE_SPARSE_INPUT"
assert status["orfs"]["candidate_orfs"]>=1
assert status["orfs"]["primary_orfs_fabricated"] is False

with (root/contract["paths"]["orfs"]).open(newline="") as handle:
    orfs=list(csv.DictReader(handle,delimiter="\t"))
assert any(row["orf_status"]=="ORF_CANDIDATE_AVAILABLE" for row in orfs)
assert not any(row["primary_orf"].lower()=="true" for row in orfs)
assert not any("non-coding" in value.lower() for row in orfs for value in row.values())

with (root/contract["paths"]["integrated"]).open(newline="") as handle:
    integrated=list(csv.DictReader(handle,delimiter="\t"))
assert {row["sample"] for row in integrated}=={"ont_real","pacbio_real"}
assert all(row["primary_orf_id"]=="NA" for row in integrated)
assert any(row["orf_status"]=="ORF_CANDIDATE_AVAILABLE" for row in integrated)

versions=contract["runtime_versions"]
assert versions["differential_expression"]["edgeR"]=="4.10.5"
assert versions["transcript_usage"]["satuRn"]=="1.20.0"
assert versions["isoform_switches"]["IsoformSwitchAnalyzeR"]=="2.12.0"
assert versions["enrichment"]["clusterProfiler"]=="4.20.2"
assert versions["SUPPA"]["release"]=="2.4" and versions["SUPPA"]["cli_reported"]=="2.3"
assert versions["TransDecoder"]["TransDecoder"]=="6.0.0"
print("Block 4 real compact-fixture contract passed.")
