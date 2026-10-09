#!/usr/bin/env python3
"""Validate the normalized Block 1 record and its reference identity."""

import argparse
import csv
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for chunk in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(chunk)
    return digest.hexdigest()


def read_validation(path: Path) -> dict[str, str]:
    with path.open(newline="") as handle:
        rows = list(csv.reader(handle, delimiter="\t"))
    if not rows or rows[0] != ["check", "value"]:
        raise ValueError("reference validation report has an invalid header")
    return {row[0]: row[1] for row in rows[1:] if len(row) >= 2}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--reference-validation", required=True, type=Path)
    parser.add_argument("--fasta", required=True, type=Path)
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--library-type", required=True)
    parser.add_argument("--genome-build", required=True)
    parser.add_argument("--annotation-release", required=True)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()

    manifest = json.loads(args.manifest.read_text())
    required = {
        "contract_version",
        "sample",
        "platform",
        "library_type",
        "condition",
        "replicate",
        "original_input",
        "sorted_bam",
        "bam_index",
        "alignment_qc_directory",
        "reference",
        "preprocessing",
        "alignment",
        "runtime_tool_versions",
        "reference_validation",
    }
    missing = sorted(required.difference(manifest))
    if missing:
        raise ValueError(f"normalized manifest is missing required fields: {', '.join(missing)}")

    expected = {
        "sample": args.sample,
        "platform": args.platform,
        "library_type": args.library_type,
    }
    for field, value in expected.items():
        if manifest.get(field) != value:
            raise ValueError(
                f"normalized manifest {field} mismatch: expected {value!r}, "
                f"observed {manifest.get(field)!r}"
            )

    reference = manifest.get("reference")
    if not isinstance(reference, dict):
        raise ValueError("normalized manifest reference field must be an object")
    reference_expected = {
        "genome_build": args.genome_build,
        "annotation_release": args.annotation_release,
    }
    for field, value in reference_expected.items():
        if reference.get(field) != value:
            raise ValueError(
                f"normalized manifest reference {field} mismatch: expected {value!r}, "
                f"observed {reference.get(field)!r}"
            )

    validation = read_validation(args.reference_validation)
    if validation.get("status") != "PASS":
        raise ValueError("Block 1 reference validation did not report PASS")

    observed_fasta = sha256(args.fasta)
    observed_gtf = sha256(args.gtf)
    if validation.get("fasta_sha256") != observed_fasta:
        raise ValueError("reference FASTA checksum differs from the Block 1 validated FASTA")
    if validation.get("annotation_sha256") != observed_gtf:
        raise ValueError("annotation checksum differs from the Block 1 validated GTF")

    report = {
        "status": "PASS",
        "sample": args.sample,
        "block1_contract_version": manifest["contract_version"],
        "genome_build": args.genome_build,
        "annotation_release": args.annotation_release,
        "fasta_sha256": observed_fasta,
        "annotation_sha256": observed_gtf,
    }
    args.output.write_text(json.dumps(report, indent=2) + "\n")


if __name__ == "__main__":
    main()
