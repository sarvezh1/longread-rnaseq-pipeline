#!/usr/bin/env python3
"""Extract one production IsoQuant sample into a stable Block 4 bundle."""

import argparse
import csv
import json
import shutil
from pathlib import Path


def read_values(path: Path, value_name: str) -> dict[str, float]:
    if not path.exists():
        return {}
    with path.open(newline="") as handle:
        rows = csv.DictReader(handle, delimiter="\t")
        result = {}
        for row in rows:
            feature = row.get("feature_id", "")
            if feature and not feature.startswith("__"):
                result[feature] = float(row[value_name])
        return result


def write_values(path: Path, values: dict[str, float], field: str, integer: bool = False) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature_id", field])
        for feature in sorted(values):
            value = values[feature]
            if integer and abs(value - round(value)) > 1e-8:
                raise ValueError(f"Raw count for {feature} is not integer-valued: {value}")
            writer.writerow([feature, int(round(value)) if integer else f"{value:.10g}"])


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--library-type", required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--isoquant-dir", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--structure", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--contract", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    sample_dir = args.isoquant_dir / args.sample
    prefix = sample_dir / args.sample
    output = args.output
    output.mkdir(parents=True, exist_ok=True)

    gene_counts = read_values(Path(f"{prefix}.discovered_gene_counts.tsv"), "count")
    transcript_counts = read_values(Path(f"{prefix}.discovered_transcript_counts.tsv"), "count")
    gene_tpm = read_values(Path(f"{prefix}.discovered_gene_tpm.tsv"), "TPM")
    transcript_tpm = read_values(Path(f"{prefix}.discovered_transcript_tpm.tsv"), "TPM")
    # Legitimate empty IsoQuant outputs may omit TPM files entirely.
    for feature in gene_counts:
        gene_tpm.setdefault(feature, 0.0)
    for feature in transcript_counts:
        transcript_tpm.setdefault(feature, 0.0)

    write_values(output / "gene_counts.tsv", gene_counts, "count", integer=True)
    write_values(output / "transcript_counts.tsv", transcript_counts, "count", integer=True)
    write_values(output / "gene_tpm.tsv", gene_tpm, "TPM")
    write_values(output / "transcript_tpm.tsv", transcript_tpm, "TPM")

    for source, name in [
        (args.gtf, "transcripts.gtf"), (args.fasta, "transcripts.fasta"),
        (args.structure, "transcript_structure.tsv"),
        (args.classification, "sqanti3_classification.tsv"),
        (args.contract, "block2_contract.json"),
    ]:
        shutil.copyfile(source, output / name)

    metadata = {
        "sample": args.sample, "condition": args.condition, "replicate": args.replicate,
        "platform": args.platform, "library_type": args.library_type,
        "production_caller": "IsoQuant", "block2_contract": args.contract.name,
        "gene_count_total": sum(gene_counts.values()),
        "transcript_count_total": sum(transcript_counts.values()),
    }
    (output / "sample_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


if __name__ == "__main__":
    main()
