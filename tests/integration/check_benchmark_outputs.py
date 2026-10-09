#!/usr/bin/env python3
"""Independent integrity checks for real Block 3 compact-fixture outputs."""

import csv
import json
import sys
from pathlib import Path


def rows(path: Path) -> list[dict]:
    assert path.is_file(), path
    with path.open(newline="") as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def fasta_ids(path: Path) -> list[str]:
    assert path.is_file(), path
    return [line[1:].split()[0] for line in path.read_text().splitlines() if line.startswith(">")]


def main() -> None:
    root = Path(sys.argv[1])
    expected = {"ont_real": "ont", "pacbio_real": "pacbio"}
    callers = ("isoquant", "flair", "bambu")
    for sample, platform in expected.items():
        summary_path = root / "benchmark" / "summary" / f"{sample}.benchmark_summary.json"
        contract_path = root / "benchmark" / "summary" / f"{sample}.benchmark_contract.json"
        summary = json.loads(summary_path.read_text())
        contract = json.loads(contract_path.read_text())
        assert summary["sample"] == sample and summary["platform"] == platform
        assert contract["sample"] == sample and contract["platform"] == platform
        assert contract["default_production_caller"] == "isoquant"
        assert contract["common_sqanti3_characterization"] is True
        assert contract["caller_ids_are_identity"] is False
        assert contract["consensus_transcriptome_created"] is False
        assert set(contract["runtime_provenance"]) == {
            f"{sample}.isoquant.versions.yml", f"{sample}.flair.versions.yml", f"{sample}.bambu.versions.yml",
            f"{sample}.isoquant.sqanti3.versions.yml", f"{sample}.flair.sqanti3.versions.yml", f"{sample}.bambu.sqanti3.versions.yml",
        }
        for caller in callers:
            structures = root / "benchmark" / "concordance" / "transcripts" / f"{sample}.{caller}.standardized_transcripts.tsv"
            records = rows(structures)
            assert len(records) == summary["caller_counts"][caller]["transcripts"]
            assert all(row["sample"] == sample and row["platform"] == platform and row["caller"] == caller for row in records)
            assert all(row["caller_transcript_id"] for row in records)
            sqanti = root / "benchmark" / "sqanti3" / caller / f"{sample}.{caller}.sqanti3"
            classification = sqanti / f"{sample}.{caller}.sqanti3_classification.txt"
            transcript_fasta = sqanti / f"{sample}.{caller}.sqanti3_corrected.fasta"
            class_rows = rows(classification)
            assert sorted(row["caller_transcript_id"] for row in records) == sorted(row["isoform"] for row in class_rows)
            assert sorted(fasta_ids(transcript_fasta)) == sorted(row["isoform"] for row in class_rows)
        membership = rows(root / "benchmark" / "concordance" / "membership" / f"{sample}.membership.tsv")
        assert sum(summary["membership_counts"].values()) == len(membership)
        pairwise = rows(root / "benchmark" / "concordance" / "intron_chains" / f"{sample}.pairwise_metrics.tsv")
        assert len(pairwise) == 12
        assert all(row["status"] in {"defined", "one_empty", "both_empty"} for row in pairwise)

    ont = json.loads((root / "benchmark" / "summary" / "ont_real.benchmark_summary.json").read_text())
    pacbio = json.loads((root / "benchmark" / "summary" / "pacbio_real.benchmark_summary.json").read_text())
    assert all(ont["caller_counts"][caller]["transcripts"] == 0 for caller in callers)
    assert pacbio["caller_counts"]["isoquant"]["transcripts"] == 1
    assert pacbio["caller_counts"]["flair"]["transcripts"] == 0
    assert pacbio["caller_counts"]["bambu"]["transcripts"] == 0
    assert pacbio["membership_counts"] == {"isoquant_only": 1}
    for sample in expected:
        flair_version = (root / "provenance" / "software_versions" / f"{sample}.flair.versions.yml").read_text()
        bambu_version = (root / "provenance" / "software_versions" / f"{sample}.bambu.versions.yml").read_text()
        assert "flair: 3.0.1" in flair_version and "NO_USABLE_SPLICE_JUNCTIONS" in flair_version
        assert "bambu: 3.14.0" in bambu_version and "bioconductor: 3.23" in bambu_version and "r: 4.6.1" in bambu_version
    print("Block 3 real benchmark outputs passed independent integrity checks.")


if __name__ == "__main__":
    main()
