#!/usr/bin/env python3
"""Validate the frozen v1 sample-sheet contract without running Nextflow."""

import argparse
import csv
from pathlib import Path

COLUMNS = ["sample", "platform", "library_type", "preprocessing", "condition", "replicate", "reads", "primer_fasta", "primer_pair", "ont_kit", "require_polya"]
KITS = {"PCS109", "PCS110", "PCS111", "PCS114", "PCB111", "PCB114"}


def boolean(value: str, field: str, sample: str) -> bool:
    if value.strip().lower() not in {"true", "false"}:
        raise ValueError(f"Sample '{sample}' has invalid {field} '{value}'; use true or false")
    return value.strip().lower() == "true"


def validate(path: Path, allow_missing_files: bool = False) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8-sig") as handle:
        reader = csv.DictReader(handle)
        if reader.fieldnames != COLUMNS:
            raise ValueError(f"Expected columns in this exact order: {', '.join(COLUMNS)}")
        rows = list(reader)
    if not rows:
        raise ValueError("Sample sheet contains no records")
    seen = set()
    for row in rows:
        sample = row["sample"].strip()
        for field in ("sample", "platform", "library_type", "preprocessing", "condition", "replicate", "reads"):
            if not row[field].strip():
                raise ValueError(f"Sample '{sample or '<unknown>'}' has empty required field '{field}'")
        if sample in seen:
            raise ValueError(f"Duplicate sample identifier: {sample}")
        seen.add(sample)
        if not sample[0].isalnum() or any(char not in "abcdefghijklmnopqrstuvwxyzABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789_.-" for char in sample):
            raise ValueError(f"Invalid sample identifier: {sample}")
        platform, library = row["platform"].strip(), row["library_type"].strip()
        preprocessing = row["preprocessing"].strip()
        require_polya = boolean(row["require_polya"], "require_polya", sample)
        reads = row["reads"].strip()
        if platform == "ont":
            if library not in {"cdna", "pcr_cdna"} or preprocessing not in {"pychopper", "none"}:
                raise ValueError(f"Unsupported ONT configuration for sample '{sample}'")
            if not reads.lower().endswith((".fastq", ".fq", ".fastq.gz", ".fq.gz")):
                raise ValueError(f"ONT sample '{sample}' requires FASTQ input")
            if preprocessing == "pychopper" and row["ont_kit"].strip() not in KITS:
                raise ValueError(f"ONT Pychopper sample '{sample}' requires a supported ont_kit")
            if row["primer_fasta"].strip() or row["primer_pair"].strip() or require_polya:
                raise ValueError(f"ONT sample '{sample}' sets a PacBio-only field")
        elif platform == "pacbio":
            if not reads.lower().endswith(".bam"):
                raise ValueError(f"PacBio sample '{sample}' requires BAM input")
            if row["ont_kit"].strip():
                raise ValueError(f"PacBio sample '{sample}' must not set ont_kit")
            if library == "isoseq_hifi":
                if preprocessing != "isoseq" or not row["primer_fasta"].strip():
                    raise ValueError(f"PacBio HiFi sample '{sample}' requires isoseq preprocessing and primer_fasta")
            elif library == "isoseq_flnc":
                if preprocessing != "none" or row["primer_fasta"].strip() or row["primer_pair"].strip() or require_polya:
                    raise ValueError(f"PacBio FLNC sample '{sample}' must bypass lima/refine fields")
            else:
                raise ValueError(f"Unsupported PacBio library type for sample '{sample}'")
        else:
            raise ValueError(f"Unsupported platform for sample '{sample}': {platform}")
        if not allow_missing_files:
            for field in ("reads", "primer_fasta"):
                value = row[field].strip()
                if value and not Path(value).exists():
                    raise ValueError(f"Sample '{sample}' {field} does not exist: {value}")
    return rows


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("sample_sheet", type=Path)
    parser.add_argument("--allow-missing-files", action="store_true")
    args = parser.parse_args()
    rows = validate(args.sample_sheet, args.allow_missing_files)
    print(f"Validated {len(rows)} sample records: {args.sample_sheet}")


if __name__ == "__main__":
    main()

