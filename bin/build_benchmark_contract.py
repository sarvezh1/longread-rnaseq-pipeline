#!/usr/bin/env python3
"""Build the machine-readable Block 3 benchmark contract and runtime record."""

import argparse
import hashlib
import json
from pathlib import Path


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--summary", required=True, type=Path)
    parser.add_argument("--manifest", required=True, type=Path)
    parser.add_argument("--version", action="append", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    summary = json.loads(args.summary.read_text())
    if summary.get("sample") != args.sample or summary.get("platform") != args.platform:
        raise ValueError("benchmark summary metadata does not match the contract request")
    runtime = {path.name: {"sha256": sha256(path), "observed": path.read_text().splitlines()} for path in args.version}
    contract = {
        "contract_version": "3.0",
        "sample": args.sample,
        "platform": args.platform,
        "condition": args.condition,
        "replicate": args.replicate,
        "input_contract": {"block1_manifest": args.manifest.name, "sha256": sha256(args.manifest)},
        "callers": ["isoquant", "flair", "bambu"],
        "default_production_caller": "isoquant",
        "common_sqanti3_characterization": True,
        "caller_ids_are_identity": False,
        "consensus_transcriptome_created": False,
        "comparison_summary": summary,
        "published_paths": {
            "caller_outputs": {caller: f"benchmark/{caller}" for caller in ("isoquant", "flair", "bambu")},
            "sqanti3": {caller: f"benchmark/sqanti3/{caller}" for caller in ("isoquant", "flair", "bambu")},
            "transcripts": "benchmark/concordance/transcripts",
            "intron_chains": "benchmark/concordance/intron_chains",
            "splice_junctions": "benchmark/concordance/splice_junctions",
            "transcript_ends": "benchmark/concordance/transcript_ends",
            "membership": "benchmark/concordance/membership",
            "summary": "benchmark/summary",
        },
        "runtime_provenance": runtime,
        "future_hooks": ["replicate_reproducibility", "platform_comparison", "caller_by_platform_interaction", "explicit_consensus_policy"],
    }
    args.output.write_text(json.dumps(contract, indent=2) + "\n")


if __name__ == "__main__":
    main()
