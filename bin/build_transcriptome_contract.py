#!/usr/bin/env python3
"""Create the standardized per-sample Block 2 transcriptome manifest."""

import argparse
import json
from pathlib import Path


def main() -> None:
    parser = argparse.ArgumentParser()
    for name in (
        "sample", "platform", "library_type", "condition", "replicate",
        "genome_build", "annotation_release", "fasta_sha256", "annotation_sha256",
        "isoquant_data_type", "isoquant_full_length", "isoquant_novel_unspliced",
        "primary_gtf", "transcript_fasta", "classification", "junctions",
        "support_table", "counts_directory", "sqanti_directory", "structure_table",
        "summary_json", "block1_manifest", "isoquant_versions", "sqanti3_versions",
        "output",
    ):
        parser.add_argument(f"--{name.replace('_', '-')}", required=True)
    args = parser.parse_args()

    manifest = {
        "contract_version": "2.0",
        "sample": args.sample,
        "platform": args.platform,
        "library_type": args.library_type,
        "condition": args.condition,
        "replicate": args.replicate,
        "reference": {
            "genome_build": args.genome_build,
            "annotation_release": args.annotation_release,
            "fasta_sha256": args.fasta_sha256,
            "annotation_sha256": args.annotation_sha256,
        },
        "primary_discovery": {
            "engine": "IsoQuant",
            "data_type": args.isoquant_data_type,
            "bulk_mode": True,
            "full_length_input": args.isoquant_full_length.lower() == "true",
            "novel_unspliced_policy": args.isoquant_novel_unspliced,
            "filtering": "none",
        },
        "structural_qc": {
            "engine": "SQANTI3",
            "classification_schema": "SQANTI3 6.x",
            "orf_prediction": "skipped",
        },
        "outputs": {
            "primary_transcript_annotation": args.primary_gtf,
            "transcript_fasta": args.transcript_fasta,
            "transcript_classification": args.classification,
            "splice_junction_classification": args.junctions,
            "transcript_support": args.support_table,
            "isoquant_output_directory": args.counts_directory,
            "sqanti3_output_directory": args.sqanti_directory,
            "structural_identity_table": args.structure_table,
            "structural_summary": args.summary_json,
        },
        "provenance": {
            "block1_normalized_manifest": args.block1_manifest,
            "isoquant_versions": args.isoquant_versions,
            "sqanti3_versions": args.sqanti3_versions,
        },
        "identity_note": (
            "Transcript identity is defined by chromosome, strand, exon coordinates, "
            "intron chain, transcript start, and transcript end; transcript IDs are run-local labels."
        ),
        "curated_transcriptome": None,
    }
    Path(args.output).write_text(json.dumps(manifest, indent=2) + "\n")


if __name__ == "__main__":
    main()
