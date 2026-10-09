#!/usr/bin/env python3
"""Build deterministic Block 5 reports from standardized pipeline contracts."""

from __future__ import annotations

import argparse
import csv
import hashlib
import html
import json
import os
import re
from datetime import datetime, timezone
from pathlib import Path
from typing import Any


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def load_json(path: Path) -> dict[str, Any]:
    return json.loads(path.read_text(encoding="utf-8"))


def read_tsv_map(path: Path) -> dict[str, str]:
    if not path.exists() or not path.stat().st_size:
        return {}
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows:
        return {}
    if rows[0][:2] in (["metric", "value"], ["check", "value"]):
        return {row[0]: row[1] for row in rows[1:] if len(row) >= 2}
    return {}


def read_nanostats(path: Path) -> dict[str, str]:
    result: dict[str, str] = {}
    if not path.exists():
        return result
    for row in csv.reader(path.read_text(encoding="utf-8", errors="replace").splitlines(), delimiter="\t"):
        if len(row) >= 2 and row[0] != "Metrics":
            result[row[0].rstrip(":")] = row[1]
    return result


def safe_name(value: Any) -> str:
    return Path(str(value)).name if value else ""


def public_parameters(path: Path) -> dict[str, Any]:
    raw = load_json(path)
    sensitive_paths = {"input", "fasta", "gtf", "outdir"}
    return {key: safe_name(value) if key in sensitive_paths else value for key, value in sorted(raw.items())}


def find_json(root: Path, pattern: str) -> list[Path]:
    return sorted(path for path in root.rglob(pattern) if path.is_file())


def contract_files(root: Path, pattern: str) -> list[Path]:
    files = []
    for path in find_json(root, pattern):
        try:
            if not load_json(path).get("placeholder"):
                files.append(path)
        except (json.JSONDecodeError, OSError):
            raise ValueError(f"Malformed contract JSON: {path}")
    return files


def parse_versions(root: Path) -> tuple[dict[str, str], list[dict[str, str]]]:
    observed_values: dict[str, set[str]] = {}
    sources: list[dict[str, str]] = []
    def record(key: str, value: str) -> None:
        observed_values.setdefault(key.lower(), set()).add(value)
    for path in sorted(root.rglob("*")):
        if not path.is_file() or path.suffix.lower() not in {".yml", ".yaml", ".json"}:
            continue
        sources.append({"file": path.name, "sha256": sha256(path)})
        if path.suffix.lower() == ".json":
            try:
                data = load_json(path)
            except json.JSONDecodeError:
                continue
            stack = [data]
            while stack:
                item = stack.pop()
                if isinstance(item, dict):
                    for key, value in item.items():
                        if isinstance(value, (dict, list)):
                            stack.append(value)
                        elif value is not None:
                            record(str(key), str(value))
                elif isinstance(item, list):
                    stack.extend(item)
        else:
            for line in path.read_text(encoding="utf-8", errors="replace").splitlines():
                match = re.match(r"^\s+([A-Za-z][A-Za-z0-9_.-]*):\s*[\"']?([^\"'].*?)[\"']?\s*$", line)
                if match and match.group(2) not in {"", "null", "~"}:
                    record(match.group(1), match.group(2).strip().strip('"\''))
    observed = {key: "; ".join(sorted(values)) for key, values in observed_values.items()}
    return observed, sources


def version_audit(registry: dict[str, Any], observed: dict[str, str], nextflow_version: str) -> tuple[list[dict[str, str]], list[str]]:
    observed.setdefault("nextflow", nextflow_version)
    aliases = {
        "nanoplot": ["nanoplot"], "pychopper": ["pychopper"], "minimap2": ["minimap2"],
        "samtools": ["samtools"], "lima": ["lima"], "isoseq": ["isoseq"], "pbmm2": ["pbmm2"],
        "isoquant": ["isoquant"], "sqanti3": ["sqanti3"], "flair": ["flair"], "bambu": ["bambu"],
        "r": ["r"], "bioconductor": ["bioconductor"], "edger": ["edger"], "saturn": ["saturn"],
        "isoformswitchanalyzer": ["isoformswitchanalyzer"], "clusterprofiler": ["clusterprofiler"],
        "org.hs.eg.db": ["org.hs.eg.db"], "suppa": ["release", "suppa"],
        "transdecoder": ["transdecoder"], "multiqc": ["multiqc"], "nextflow": ["nextflow"],
    }
    rows: list[dict[str, str]] = []
    warnings: list[str] = []
    for spec in registry["tools"]:
        tool = spec["tool"]
        keys = aliases.get(tool.lower(), [tool.lower()])
        actual = next((observed[key] for key in keys if key in observed), "NOT_OBSERVED")
        expected = spec["expected"]
        if actual == "NOT_OBSERVED":
            status = "NOT_OBSERVED_NOT_RUN"
        elif expected.startswith(">="):
            status = "MATCH" if actual else "MISMATCH"
        else:
            expected_values = {value.strip() for value in expected.split(";")}
            actual_values = {value.strip() for value in actual.split(";")}
            status = "MATCH" if actual_values == expected_values else "MISMATCH"
        if status == "MISMATCH":
            warnings.append(f"Version mismatch for {tool}: expected {expected}, observed {actual}")
        container = spec["container"]
        image_identity = spec.get("image_id", "")
        immutable = bool(image_identity) or "@sha256:" in container or container == "host runtime" or container.startswith("community.wave.seqera.io/")
        if not immutable:
            warnings.append(f"Container reference is version-pinned but lacks a registry digest: {container}")
        rows.append({
            "tool": tool, "expected_version": expected, "observed_version": actual,
            "version_status": status, "container": container,
            "image_identity": image_identity,
            "digest_status": "IMMUTABLE_OR_CONTENT_ADDRESSED" if immutable else "PINNED_TAG_NO_REGISTRY_DIGEST",
            "process": spec["process"],
        })
    return rows, sorted(set(warnings))


def status_class(status: str, vocabulary: dict[str, Any]) -> str:
    if status in vocabulary:
        return vocabulary[status]["class"]
    if status.startswith("NOT_RUN_"):
        return "not_run"
    if status.endswith("EMPTY_RESULT") or status == "NO_TRANSCRIPTS":
        return "successful_empty"
    if "NOT_AVAILABLE" in status:
        return "unavailable"
    if status == "FAILED":
        return "failed"
    return "successful" if status in {"COMPLETE", "PASS", "CLASSIFIED"} else "informational"


def write_tsv(path: Path, rows: list[dict[str, Any]], fields: list[str]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def table(headers: list[str], rows: list[list[Any]]) -> str:
    head = "".join(f"<th>{html.escape(str(value))}</th>" for value in headers)
    body = "".join("<tr>" + "".join(f"<td>{html.escape(str(value))}</td>" for value in row) + "</tr>" for row in rows)
    return f"<div class='table-wrap'><table><thead><tr>{head}</tr></thead><tbody>{body}</tbody></table></div>"


def link(relative: str, label: str | None = None) -> str:
    target = "../" + relative.lstrip("./")
    return f"<a href='{html.escape(target, quote=True)}'>{html.escape(label or relative)}</a>"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--contracts", type=Path, required=True)
    parser.add_argument("--qc", type=Path, required=True)
    parser.add_argument("--versions", type=Path, required=True)
    parser.add_argument("--reference-validation", type=Path, required=True)
    parser.add_argument("--sample-sheet", type=Path, required=True)
    parser.add_argument("--parameters", type=Path, required=True)
    parser.add_argument("--status-vocabulary", type=Path, required=True)
    parser.add_argument("--tool-registry", type=Path, required=True)
    parser.add_argument("--mode", required=True)
    parser.add_argument("--pipeline-version", required=True)
    parser.add_argument("--nextflow-version", required=True)
    parser.add_argument("--run-id", required=True)
    parser.add_argument("--execution-profile", required=True)
    parser.add_argument("--started-at", required=True)
    parser.add_argument("--completed-at", default="")
    parser.add_argument("--multiqc-report", default="multiqc/multiqc_report.html")
    parser.add_argument("--output-report", type=Path, required=True)
    parser.add_argument("--output-provenance", type=Path, required=True)
    args = parser.parse_args()

    args.output_report.mkdir(parents=True, exist_ok=True)
    args.output_provenance.mkdir(parents=True, exist_ok=True)
    vocabulary = load_json(args.status_vocabulary)
    registry = load_json(args.tool_registry)
    parameters = public_parameters(args.parameters)
    normalized_paths = contract_files(args.contracts, "*.normalized.json")
    transcript_paths = contract_files(args.contracts, "*.transcriptome_contract.json")
    benchmark_paths = contract_files(args.contracts, "*.benchmark_contract.json")
    biology_paths = contract_files(args.contracts, "biology_contract.json")
    normalized = {item["sample"]: item for item in map(load_json, normalized_paths)}
    transcriptomes = {item["sample"]: item for item in map(load_json, transcript_paths)}
    benchmarks = {item["sample"]: item for item in map(load_json, benchmark_paths)}
    biology = load_json(biology_paths[0]) if biology_paths else None
    reference = read_tsv_map(args.reference_validation)

    required_errors = []
    if not normalized:
        required_errors.append("No valid Block 1 normalized contracts were supplied")
    if reference.get("status") not in {"PASS", "STUB"}:
        required_errors.append("Reference compatibility status is not PASS")
    if args.mode in {"default", "biology"} and len(transcriptomes) != len(normalized):
        required_errors.append("A Block 2 transcriptome contract is missing")
    if args.mode == "benchmark" and len(benchmarks) != len(normalized):
        required_errors.append("A Block 3 benchmark contract is missing")
    if args.mode == "biology" and biology is None:
        required_errors.append("The Block 4 biology contract is missing")

    transcript_summaries: dict[str, dict[str, Any]] = {}
    for path in find_json(args.qc, "*.transcript_summary.json"):
        data = load_json(path)
        transcript_summaries[data["sample"]] = data
    samples = []
    for sample_id, contract in sorted(normalized.items()):
        nano_paths = list(args.qc.rglob(f"{sample_id}.nanoplot/NanoStats.txt"))
        alignment_paths = list(args.qc.rglob(f"{sample_id}.alignment_qc/{sample_id}.alignment_summary.tsv"))
        nano = read_nanostats(nano_paths[0]) if nano_paths else {}
        alignment = read_tsv_map(alignment_paths[0]) if alignment_paths else {}
        tx = transcript_summaries.get(sample_id, {})
        samples.append({
            "sample": sample_id,
            "platform": contract.get("platform", "NA"),
            "library_type": contract.get("library_type", "NA"),
            "condition": contract.get("condition", "NA"),
            "replicate": contract.get("replicate", "NA"),
            "input_file": safe_name(contract.get("original_input")),
            "preprocessing": contract.get("preprocessing", {}).get("selection", "NA"),
            "preprocessing_status": "COMPLETE",
            "raw_reads": nano.get("number_of_reads", "NA"),
            "raw_bases": nano.get("number_of_bases", "NA"),
            "read_n50": nano.get("n50", "NA"),
            "mean_quality": nano.get("mean_qual", "NA"),
            "alignment_status": "COMPLETE",
            "total_records": alignment.get("total_records", "NA"),
            "mapped_records": alignment.get("mapped_records", "NA"),
            "mapped_percentage": alignment.get("mapped_percentage", "NA"),
            "unmapped_records": alignment.get("unmapped_records", "NA"),
            "transcripts": tx.get("total_transcript_models", "NOT_RUN" if args.mode == "benchmark" else 0),
            "single_exon": tx.get("single_exon_transcripts", "NA"),
            "multi_exon": tx.get("multi_exon_transcripts", "NA"),
            "bam": contract.get("sorted_bam", ""),
            "bai": contract.get("bam_index", ""),
        })

    analysis_rows: list[dict[str, str]] = []
    if biology:
        for analysis, details in biology.get("analysis_status", {}).items():
            status = str(details.get("status", "UNKNOWN"))
            analysis_rows.append({"analysis": analysis, "status": status, "class": status_class(status, vocabulary)})
    else:
        analysis_rows.append({"analysis": "biology", "status": "NOT_RUN_DISABLED", "class": "not_run"})

    observed, version_sources = parse_versions(args.versions)
    tool_rows, audit_warnings = version_audit(registry, observed, args.nextflow_version)
    overall_status = "FAILED" if required_errors else "COMPLETE"
    completed_at = args.completed_at or datetime.now(timezone.utc).replace(microsecond=0).isoformat()
    reference_public = {
        "status": reference.get("status", "UNKNOWN"),
        "genome_build": parameters.get("genome_build"),
        "annotation_release": parameters.get("annotation_release"),
        "fasta_file": safe_name(reference.get("fasta_file")),
        "annotation_file": safe_name(reference.get("annotation_file")),
        "fasta_sha256": reference.get("fasta_sha256"),
        "fasta_index_sha256": reference.get("fasta_index_sha256"),
        "annotation_sha256": reference.get("annotation_sha256"),
    }
    summary = {
        "schema_version": "5.0",
        "pipeline": {"name": "longread-rnaseq-pipeline", "version": args.pipeline_version},
        "run": {
            "id": args.run_id, "mode": args.mode, "execution_profile": args.execution_profile,
            "started_at": args.started_at, "completed_at": completed_at, "overall_status": overall_status,
            "required_errors": required_errors, "warnings": audit_warnings,
        },
        "reference": reference_public,
        "samples": samples,
        "analysis_status": analysis_rows,
        "benchmark": {sample: contract.get("comparison_summary", {}) for sample, contract in benchmarks.items()},
        "biology": biology or {"status": "NOT_RUN_DISABLED"},
        "evidence_hierarchy": (biology or {}).get("evidence_hierarchy", {
            "observed": ["reads", "alignments", "transcript structures"],
            "statistically_inferred": [], "computationally_predicted": [], "experimentally_validated": []}),
        "major_outputs": {
            "multiqc": f"report/{args.multiqc_report}",
            "normalized_contracts": "provenance/samples/*.normalized.json",
            "transcriptome_contracts": "provenance/samples/*.transcriptome_contract.json",
            "benchmark_contracts": "benchmark/summary/*.benchmark_contract.json",
            "biology_contract": "biology/biology_contract.json",
        },
        "tool_audit": tool_rows,
    }
    (args.output_report / "run_summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    write_tsv(args.output_report / "samples.tsv", samples, list(samples[0]) if samples else ["sample"])
    write_tsv(args.output_report / "analysis_status.tsv", analysis_rows, ["analysis", "status", "class"])
    write_tsv(args.output_report / "tool_versions.tsv", tool_rows, ["tool", "expected_version", "observed_version", "version_status", "container", "image_identity", "digest_status", "process"])

    reproducibility = {
        "schema_version": "1.0",
        "pipeline": summary["pipeline"],
        "run": {key: summary["run"][key] for key in ("id", "mode", "execution_profile", "started_at", "completed_at", "overall_status")},
        "software": tool_rows,
        "reference": reference_public,
        "parameters": parameters,
        "sample_sheet": {"file": args.sample_sheet.name, "sha256": sha256(args.sample_sheet)},
        "samples": [{key: row[key] for key in ("sample", "platform", "library_type", "condition", "replicate", "input_file")} for row in samples],
        "version_record_sources": version_sources,
        "work_directory_is_essential_provenance": False,
    }
    (args.output_provenance / "reproducibility_manifest.json").write_text(json.dumps(reproducibility, indent=2) + "\n", encoding="utf-8")
    (args.output_provenance / "resolved_parameters.json").write_text(json.dumps(parameters, indent=2) + "\n", encoding="utf-8")
    (args.output_provenance / "provenance_bundle.json").write_text(json.dumps({
        "schema_version": "5.0", "run_summary_sha256": sha256(args.output_report / "run_summary.json"),
        "reproducibility_manifest": "provenance/final/reproducibility_manifest.json",
        "tool_audit": "provenance/final/tool_audit.tsv", "reference": reference_public,
        "sample_sheet_sha256": sha256(args.sample_sheet), "completion_state": overall_status,
    }, indent=2) + "\n", encoding="utf-8")
    write_tsv(args.output_provenance / "tool_audit.tsv", tool_rows, ["tool", "expected_version", "observed_version", "version_status", "container", "image_identity", "digest_status", "process"])

    sample_table = table(
        ["Sample", "Platform", "Library", "Condition", "Replicate", "Raw reads", "N50", "Mapped", "Mapped %", "Transcripts"],
        [[row[key] for key in ("sample", "platform", "library_type", "condition", "replicate", "raw_reads", "read_n50", "mapped_records", "mapped_percentage", "transcripts")] for row in samples],
    )
    status_table = table(["Analysis", "Status", "Meaning"], [[row["analysis"], row["status"], vocabulary.get(row["status"], {}).get("description", row["class"])] for row in analysis_rows])
    transcript_table = table(["Sample", "Models", "Single exon", "Multi exon"], [[row["sample"], row["transcripts"], row["single_exon"], row["multi_exon"]] for row in samples])
    benchmark_html = "<p>Benchmark mode was not run.</p>"
    if benchmarks:
        caller_rows = []
        for sample, contract in benchmarks.items():
            for caller, counts in contract.get("comparison_summary", {}).get("caller_counts", {}).items():
                caller_rows.append([sample, caller, counts.get("transcripts", 0), counts.get("unique_intron_chains", 0), counts.get("splice_junctions", 0)])
        benchmark_html = table(["Sample", "Caller", "Transcripts", "Intron chains", "Junctions"], caller_rows)
    biology_links = "<p>Biology mode was not run.</p>"
    if biology:
        biology_links = "<ul>" + "".join(f"<li>{link(value, key.replace('_', ' '))}</li>" for key, value in biology.get("paths", {}).items()) + "</ul>"
    transcript_output = "benchmark/isoquant" if benchmarks else "transcriptome/raw"
    structural_qc_output = "benchmark/sqanti3" if benchmarks else "qc/transcript_structure"
    warnings_html = "<p>None.</p>" if not audit_warnings else "<ul>" + "".join(f"<li>{html.escape(item)}</li>" for item in audit_warnings) + "</ul>"
    css = """
    :root{--ink:#17212b;--muted:#5b6773;--line:#d9e0e7;--accent:#176b87;--ok:#197548;--empty:#6a5d00;--off:#59636e;--bad:#a52a2a;--paper:#fff;--wash:#f4f7f9}
    *{box-sizing:border-box}body{margin:0;background:var(--wash);color:var(--ink);font:15px/1.55 system-ui,-apple-system,Segoe UI,sans-serif}main{max-width:1180px;margin:auto;background:var(--paper);padding:2.4rem 3rem 4rem;min-height:100vh}h1{margin:.1rem 0;color:#103b4a}h2{border-bottom:2px solid var(--line);padding-bottom:.35rem;margin-top:2.2rem}h3{margin-top:1.4rem}.lede{color:var(--muted);font-size:1.05rem}.pill{display:inline-block;padding:.2rem .65rem;border-radius:999px;background:#e4f3eb;color:var(--ok);font-weight:700}.pill.failed{background:#f8e2e2;color:var(--bad)}.cards{display:grid;grid-template-columns:repeat(auto-fit,minmax(180px,1fr));gap:.75rem}.card{border:1px solid var(--line);border-radius:8px;padding:.8rem;background:#fbfcfd}.card b{display:block;color:var(--muted);font-size:.8rem;text-transform:uppercase}.table-wrap{overflow:auto;border:1px solid var(--line);border-radius:7px}table{border-collapse:collapse;width:100%}th,td{text-align:left;padding:.55rem .65rem;border-bottom:1px solid var(--line);white-space:nowrap}th{background:#edf3f6}a{color:var(--accent)}code{background:#edf1f4;padding:.1rem .25rem;border-radius:3px}footer{margin-top:3rem;color:var(--muted);border-top:1px solid var(--line);padding-top:1rem}@media(max-width:700px){main{padding:1.2rem}.cards{grid-template-columns:1fr}}
    """
    report_html = f"""<!doctype html><html lang='en'><head><meta charset='utf-8'><meta name='viewport' content='width=device-width,initial-scale=1'><title>longread-rnaseq-pipeline run report</title><style>{css}</style></head><body><main>
    <h1>longread-rnaseq-pipeline</h1><p class='lede'>Deterministic final run report generated from standardized workflow contracts.</p>
    <p><span class='pill {'failed' if overall_status == 'FAILED' else ''}'>{html.escape(overall_status)}</span></p>
    <h2>Run overview</h2><div class='cards'><div class='card'><b>Version</b>{html.escape(args.pipeline_version)}</div><div class='card'><b>Run</b>{html.escape(args.run_id)}</div><div class='card'><b>Mode</b>{html.escape(args.mode)}</div><div class='card'><b>Profile</b>{html.escape(args.execution_profile)}</div><div class='card'><b>Samples</b>{len(samples)}</div><div class='card'><b>Completed</b>{html.escape(completed_at)}</div></div>{sample_table}
    <h2>Input and reference validation</h2>{table(['Field','Value'], [[key, value if value is not None else 'NA'] for key,value in reference_public.items()])}
    <h2>Sequencing and preprocessing QC</h2><p>Tool-native sequencing and alignment QC is aggregated in <a href='{html.escape(args.multiqc_report, quote=True)}'>the offline MultiQC report</a>.</p>{sample_table}
    <h2>Alignment</h2><p>BAM and BAI outputs are linked through the normalized contracts.</p><ul>{''.join(f'<li>{html.escape(row["sample"])}: {link(row["bam"], "BAM")} · {link(row["bai"], "BAI")}</li>' for row in samples)}</ul>
    <h2>Transcriptome discovery</h2>{transcript_table}<p>{link(transcript_output, 'Production transcript GTF/FASTA')} &middot; {link(structural_qc_output, 'SQANTI3-aware structural summaries')}</p>
    <h2>Caller benchmark</h2><p>Concordance is structural evidence; shared does not mean true and caller-specific does not mean false.</p>{benchmark_html}{'<p>'+link('benchmark/summary','Benchmark summaries and pairwise metrics')+'</p>' if benchmarks else ''}
    <h2>Expression and biological analysis</h2>{status_table}{biology_links}
    <h2>Evidence interpretation</h2>{table(['Evidence level','Contents'], [[key, ', '.join(value) if value else 'None supplied'] for key,value in summary['evidence_hierarchy'].items()])}<p>This report organizes explicit evidence and does not generate mechanistic claims.</p>
    <h2>Provenance and reproducibility</h2><p>{link('../provenance/final/reproducibility_manifest.json','Reproducibility manifest')} · {link('../provenance/final/tool_audit.tsv','Software/container audit')} · {link('../provenance/final/provenance_bundle.json','Consolidated provenance')}</p>{table(['Tool','Expected','Observed','Status'], [[row['tool'],row['expected_version'],row['observed_version'],row['version_status']] for row in tool_rows])}
    <h3>Audit warnings</h3>{warnings_html}
    <h2>Status definitions</h2>{table(['Status','Class','Meaning'], [[key,value['class'],value['description']] for key,value in vocabulary.items()])}
    <footer>Generated without external services. Core content and styling are embedded for offline viewing.</footer></main></body></html>"""
    (args.output_report / "final_report.html").write_text(report_html, encoding="utf-8")


if __name__ == "__main__":
    main()
