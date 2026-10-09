#!/usr/bin/env python3
"""Validate Block 5 real-report content, links, statuses, and privacy."""

import json
import re
import sys
from pathlib import Path

root=Path(sys.argv[1])
expected_mode=sys.argv[2]
report_dir=root/"report"
summary=json.loads((report_dir/"run_summary.json").read_text(encoding="utf-8"))
assert summary["schema_version"]=="5.0"
assert summary["pipeline"]["version"]=="0.1.0"
assert summary["run"]["mode"]==expected_mode
assert summary["run"]["overall_status"]=="COMPLETE"
assert summary["reference"]["status"]=="PASS"
assert {row["sample"] for row in summary["samples"]}=={"ont_real","pacbio_real"}
assert not any("/mnt/" in row.get("input_file","") or ":\\" in row.get("input_file","") for row in summary["samples"])
html_path=report_dir/"final_report.html"
text=html_path.read_text(encoding="utf-8")
assert "http://" not in text and "https://" not in text
assert "NOT_RUN" in text or expected_mode!="biology"
assert (report_dir/"multiqc/multiqc_report.html").is_file()
for href in re.findall(r"href='([^']+)'",text):
    if href.startswith("#"): continue
    target=(report_dir/href).resolve()
    assert target.exists(),f"broken report link: {href} -> {target}"
for relative in ("samples.tsv","analysis_status.tsv","tool_versions.tsv"):
    assert (report_dir/relative).is_file()
provenance=root/"provenance/final"
manifest=json.loads((provenance/"reproducibility_manifest.json").read_text())
bundle=json.loads((provenance/"provenance_bundle.json").read_text())
assert manifest["pipeline"]["version"]=="0.1.0"
assert manifest["work_directory_is_essential_provenance"] is False
assert bundle["completion_state"]=="COMPLETE"
tools={row["tool"]:row for row in summary["tool_audit"]}
assert tools["MultiQC"]["observed_version"]=="1.35"
assert tools["MultiQC"]["version_status"]=="MATCH"
if expected_mode=="biology":
    statuses={row["analysis"]:row["status"] for row in summary["analysis_status"]}
    assert statuses["orfs"]=="PREDICTION_NOT_AVAILABLE_SPARSE_INPUT"
    assert statuses["gene_differential_expression"]=="NOT_RUN_NO_EXPLICIT_CONTRAST"
    ont=next(row for row in summary["samples"] if row["sample"]=="ont_real")
    assert int(ont["transcripts"])==0
if expected_mode=="benchmark":
    assert summary["benchmark"]["ont_real"]["caller_counts"]["isoquant"]["transcripts"]==0
print(f"Block 5 {expected_mode} real-report contract passed.")
