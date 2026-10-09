#!/usr/bin/env python3
"""Deterministic regression cases for the Block 5 report generator."""

import csv
import json
import subprocess
import sys
import tempfile
from pathlib import Path

ROOT = Path(__file__).resolve().parents[2]


def write(path: Path, value: str) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(value, encoding="utf-8")


def run_case(name: str, mode: str, transcript_count: int, biology_statuses=None, benchmark=False, reference_status="PASS", omit_optional=False):
    with tempfile.TemporaryDirectory(prefix=f"block5-{name}-") as tmp:
        root = Path(tmp)
        contracts, qc, versions = root / "contracts", root / "qc", root / "versions"
        for path in (contracts, qc, versions): path.mkdir()
        normalized = {
            "contract_version":"1.0", "sample":"sample_1", "platform":"ont", "library_type":"cdna",
            "condition":"control", "replicate":"1", "original_input":"data/sample_1.fastq.gz",
            "sorted_bam":"alignments/ont/sample_1.sorted.bam", "bam_index":"alignments/ont/sample_1.sorted.bam.bai",
            "preprocessing":{"selection":"none"}
        }
        write(contracts / "sample_1.normalized.json", json.dumps(normalized))
        if mode in {"default", "biology"}:
            write(contracts / "sample_1.transcriptome_contract.json", json.dumps({"contract_version":"2.0","sample":"sample_1","outputs":{}}))
            write(qc / "sample_1.transcript_summary.json", json.dumps({"sample":"sample_1","total_transcript_models":transcript_count,"single_exon_transcripts":transcript_count,"multi_exon_transcripts":0}))
        if benchmark:
            counts={caller:{"transcripts":0,"unique_intron_chains":0,"splice_junctions":0} for caller in ("isoquant","flair","bambu")}
            write(contracts / "sample_1.benchmark_contract.json", json.dumps({"contract_version":"3.0","sample":"sample_1","comparison_summary":{"caller_counts":counts}}))
        if biology_statuses is not None:
            paths={} if omit_optional else {"orfs":"biology/orfs/orf_predictions.tsv"}
            write(contracts / "biology_contract.json", json.dumps({"contract_version":"4.0","paths":paths,"analysis_status":biology_statuses,"evidence_hierarchy":{"observed":["reads"],"statistically_inferred":[],"computationally_predicted":["ORFs"],"experimentally_validated":[]}}))
        write(qc / "sample_1.nanoplot" / "NanoStats.txt", "Metrics\tdataset\nnumber_of_reads\t10\nn50\t500\n")
        write(qc / "sample_1.alignment_qc" / "sample_1.alignment_summary.tsv", "metric\tvalue\ntotal_records\t10\nmapped_records\t8\nmapped_percentage\t80\nunmapped_records\t2\n")
        write(versions / "runtime.yml", 'MULTIQC:\n  MultiQC: "1.35"\n')
        write(versions / "alignment.yml", 'MINIMAP2_SPLICE:\n  samtools: "1.23.1"\n')
        write(versions / "qc.yml", 'SAMTOOLS_ALIGNMENT_QC:\n  samtools: "1.24"\n')
        write(root / "reference.tsv", f"check\tvalue\nstatus\t{reference_status}\nfasta_file\treference.fa\nfasta_sha256\tabc\nannotation_file\tannotation.gtf\nannotation_sha256\tdef\n")
        write(root / "samples.csv", "sample,platform,library_type,preprocessing,condition,replicate,reads,primer_fasta,primer_pair,ont_kit,require_polya\nsample_1,ont,cdna,none,control,1,data/sample_1.fastq.gz,,,,false\n")
        write(root / "params.json", json.dumps({"genome_build":"test","annotation_release":"test","mode":mode}))
        command=[sys.executable,str(ROOT/"bin/generate_final_report.py"),"--contracts",str(contracts),"--qc",str(qc),"--versions",str(versions),"--reference-validation",str(root/"reference.tsv"),"--sample-sheet",str(root/"samples.csv"),"--parameters",str(root/"params.json"),"--status-vocabulary",str(ROOT/"assets/reporting/status_vocabulary.json"),"--tool-registry",str(ROOT/"assets/reporting/tool_registry.json"),"--mode",mode,"--pipeline-version","0.1.0","--nextflow-version","26.04.6","--run-id",name,"--execution-profile","test","--started-at","2026-01-01T00:00:00Z","--output-report",str(root/"report"),"--output-provenance",str(root/"provenance")]
        subprocess.run(command, check=True)
        summary=json.loads((root/"report/run_summary.json").read_text())
        report=(root/"report/final_report.html").read_text()
        assert name in report and summary["samples"][0]["transcripts"] == ("NOT_RUN" if mode=="benchmark" else transcript_count)
        samtools = next(row for row in summary["tool_audit"] if row["tool"] == "samtools")
        assert samtools["observed_version"] == "1.23.1; 1.24" and samtools["version_status"] == "MATCH"
        return summary, report


def main():
    run_case("complete", "default", 2)
    run_case("zero-transcript", "default", 0)
    run_case("benchmark-empty", "benchmark", 0, benchmark=True)
    statuses={
        "orfs":{"status":"PREDICTION_NOT_AVAILABLE_SPARSE_INPUT"},
        "gene_differential_expression":{"status":"NOT_RUN_NO_EXPLICIT_CONTRAST"},
        "differential_transcript_usage":{"status":"NOT_RUN_INSUFFICIENT_REPLICATION"},
        "enrichment":{"status":"NOT_RUN_DISABLED"}
    }
    summary,report=run_case("sparse-and-not-run", "biology", 1, statuses)
    assert "PREDICTION_NOT_AVAILABLE_SPARSE_INPUT" in report
    assert "NOT_RUN_NO_EXPLICIT_CONTRAST" in report and "NOT_RUN_INSUFFICIENT_REPLICATION" in report
    _,missing_report=run_case("missing-optional", "biology", 1, statuses, omit_optional=True)
    assert "COMPLETE" in missing_report
    failed,_=run_case("failed-required", "default", 1, reference_status="FAIL")
    assert failed["run"]["overall_status"] == "FAILED"
    print("Block 5 deterministic reporting regressions passed.")


if __name__ == "__main__": main()
