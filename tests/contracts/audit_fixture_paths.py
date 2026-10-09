#!/usr/bin/env python3
"""Audit local fixtures required by routine CI and contract tests."""

from __future__ import annotations

import csv
import gzip
import re
import sys
from pathlib import Path


ROOT = Path(__file__).resolve().parents[2]
STRUCTURAL = ROOT / "tests/data/structural"
PATH_COLUMNS = ("reads", "primer_fasta")
PROTECTED_SUFFIXES = (".bam", ".cram", ".fastq", ".fq", ".fastq.gz", ".fq.gz")

# These directories are consumed as datasets or through variable filenames, so
# their members are not all visible as literal paths in the calling test script.
REQUIRED_GROUPS = {
    "tests/data/benchmark": {
        "bambu.classification.tsv", "bambu.gtf", "flair.classification.tsv",
        "flair.gtf", "isoquant.classification.tsv", "isoquant.gtf", "support.tsv",
    },
    "tests/data/biology/dataset": {
        "design_status.json", "gene_counts.tsv", "gene_tpm.tsv",
        "production_transcripts.fasta", "production_transcripts.gtf",
        "sample_metadata.tsv", "transcript_counts.tsv", "transcript_metadata.tsv",
        "transcript_presence.tsv", "transcript_tpm.tsv",
    },
    "tests/data/validation": {
        "empty.gtf", "plot.tsv", "prediction.gtf", "quant_predicted.tsv",
        "quant_truth.tsv", "reads.fastq", "scoring.tsv", "structures_r1.tsv",
        "structures_r2.tsv", "trace.tsv", "truth.gtf",
    },
}


def protected(path: str) -> bool:
    return path.endswith(PROTECTED_SUFFIXES)


def main() -> int:
    errors: list[str] = []
    required: set[str] = {
        "tests/data/structural/reference.fa",
        "tests/data/structural/annotation.gtf",
    }

    sheets = sorted(STRUCTURAL.glob("*.csv"))
    for sheet in sheets:
        required.add(sheet.relative_to(ROOT).as_posix())
        with sheet.open(newline="", encoding="utf-8-sig") as handle:
            for line_number, row in enumerate(csv.DictReader(handle), 2):
                for column in PATH_COLUMNS:
                    value = (row.get(column) or "").strip()
                    if not value:
                        continue
                    if "://" in value or Path(value).is_absolute():
                        errors.append(f"{sheet.relative_to(ROOT)}:{line_number}: {column} is not a repository-local fixture: {value}")
                        continue
                    required.add(Path(value).as_posix())

    for directory, names in REQUIRED_GROUPS.items():
        required.update(f"{directory}/{name}" for name in names)

    # Catch additional literal fixture paths introduced in CI, contract tests,
    # or the stub profile without requiring this inventory to be updated first.
    sources = list((ROOT / ".github/workflows").glob("*.yml"))
    sources += list((ROOT / "tests/contracts").glob("*.sh"))
    sources += list((ROOT / "tests/contracts").glob("*.py"))
    sources += [ROOT / "conf/test_stub.config", ROOT / "validation/nextflow.config"]
    literal = re.compile(r"tests/data/[A-Za-z0-9_./-]+")
    for source in sources:
        text = source.read_text(encoding="utf-8")
        for match in literal.findall(text):
            value = match.rstrip("./")
            if value and "." in Path(value).name:
                required.add(value)

    ignore_lines = {
        line.strip() for line in (ROOT / ".gitignore").read_text(encoding="utf-8").splitlines()
        if line.strip() and not line.lstrip().startswith("#")
    }
    for relative in sorted(required):
        path = ROOT / relative
        if not path.is_file():
            errors.append(f"Missing required checkout fixture: {relative}")
            continue
        if protected(relative) and f"!{relative}" not in ignore_lines:
            errors.append(f"Protected fixture lacks exact .gitignore exception: !{relative}")
        if relative.endswith(".bam"):
            try:
                with gzip.open(path, "rb") as handle:
                    magic = handle.read(4)
            except (OSError, EOFError) as exc:
                errors.append(f"Invalid BAM fixture {relative}: {exc}")
            else:
                if magic != b"BAM\x01":
                    errors.append(f"Invalid BAM fixture {relative}: missing BAM magic")

    if errors:
        print("Fixture-path audit failed:")
        for error in sorted(set(errors)):
            print(f"- {error}")
        return 1
    print(
        f"Fixture-path audit passed: {len(sheets)} structural sample sheets and "
        f"{len(required)} required local fixtures checked."
    )
    return 0


if __name__ == "__main__":
    sys.exit(main())
