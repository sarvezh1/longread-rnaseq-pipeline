#!/usr/bin/env python3
"""Read-only release-readiness checks for the public source tree."""

from __future__ import annotations

import json
import os
import argparse
import re
import sys
from pathlib import Path

sys.dont_write_bytecode = True
from validate_samplesheet import validate

ROOT = Path(__file__).resolve().parents[1]
SKIP_DIRS = {".git", ".nextflow", "work", "results", ".integration-fixtures", "__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache"}
REQUIRED = [
    "main.nf", "nextflow.config", "nextflow_schema.json", "README.md", "PLAN.md", "CHANGELOG.md",
    "CITATION.cff", "LICENSE", ".gitignore", "docs/architecture.md", "docs/input.md", "docs/outputs.md",
    "bin/generate_final_report.py", "assets/reporting/status_vocabulary.json", "assets/reporting/tool_registry.json",
    "docs/validation.md", "validation/workflow.nf", "validation/manifests/source_datasets.tsv",
    "validation/manifests/truth_resources.tsv", "validation/manifests/metric_definitions.json",
    "bin/validation_metrics.py", "bin/build_validation_report.py",
    ".github/workflows/ci.yml", ".github/workflows/real-integration.yml",
]
# Construct the hygiene terms without reproducing them verbatim in public source.
PROHIBITED = [
    "coding" + r" assistant",
    "generated" + r"[- ]by",
    "development" + r" conversation",
    "internal" + r" prompt",
    "chat" + "gpt",
    "code" + "x",
    "automated" + r" ai review",
]
TEXT_SUFFIXES = {".md", ".nf", ".config", ".py", ".R", ".sh", ".json", ".yaml", ".yml", ".cff", ".csv", ".txt", ".Dockerfile"}


def files():
    for directory, names, filenames in os.walk(ROOT):
        names[:] = [name for name in names if name not in SKIP_DIRS and not name.startswith(".block")]
        for name in filenames:
            path=Path(directory)/name
            try:
                if path.is_file(): yield path
            except OSError:
                continue


def local_links(path: Path, text: str):
    for target in re.findall(r"\[[^]]+\]\(([^)]+)\)", text):
        target = target.split("#", 1)[0]
        if not target or "://" in target or target.startswith(("mailto:", "#")):
            continue
        yield (path.parent / target).resolve(), target


def main() -> int:
    parser=argparse.ArgumentParser()
    parser.add_argument("--allow-git-metadata", action="store_true", help="Allow CI checkout metadata while retaining all source checks")
    args=parser.parse_args()
    errors, warnings = [], []
    for required in REQUIRED:
        if not (ROOT / required).exists(): errors.append(f"Missing required project file: {required}")
    if (ROOT / ".git").exists() and not args.allow_git_metadata: errors.append("A .git directory is present; release preparation requires a source tree without repository metadata")
    for temporary in sorted(ROOT.glob(".block*_work")) + sorted(ROOT.glob(".block*_test")) + sorted(ROOT.glob(".block*_results")):
        errors.append(f"Temporary test work directory present: {temporary.name}")
    for temporary in [ROOT / ".nextflow", ROOT / "work", ROOT / "results"]:
        if temporary.exists():
            errors.append(f"Generated workflow directory present: {temporary.name}")
    for temporary in sorted(ROOT.glob(".nextflow.log*")):
        errors.append(f"Generated workflow log present: {temporary.name}")
    for cache in ("__pycache__", ".pytest_cache", ".mypy_cache", ".ruff_cache", ".block4_test"):
        matches = [path for path in ROOT.rglob(cache) if not any(part in SKIP_DIRS for part in path.relative_to(ROOT).parts[:-1])]
        if matches: errors.append(f"Generated cache or temporary artifact present: {matches[0].relative_to(ROOT)}")
    try:
        schema=json.loads((ROOT/"nextflow_schema.json").read_text(encoding="utf-8"))
        json.loads((ROOT/"assets/reporting/status_vocabulary.json").read_text(encoding="utf-8"))
        json.loads((ROOT/"assets/reporting/tool_registry.json").read_text(encoding="utf-8"))
    except (json.JSONDecodeError, OSError) as exc:
        errors.append(f"Malformed required JSON: {exc}")
        schema={"properties":{}}
    config=(ROOT/"conf/base.config").read_text(encoding="utf-8")
    params_match=re.search(r"^params\s*\{(.*?)^\}", config, re.MULTILINE|re.DOTALL)
    active=set(re.findall(r"^\s{4}([A-Za-z][A-Za-z0-9_]*)\s*=", params_match.group(1), re.MULTILINE)) if params_match else set()
    missing_schema=sorted(active-set(schema.get("properties",{})))
    if missing_schema: errors.append("Active parameters absent from schema: "+", ".join(missing_schema))
    for sheet in sorted((ROOT/"assets/samplesheets").glob("*.csv")):
        try: validate(sheet, allow_missing_files=True)
        except ValueError as exc: errors.append(f"Invalid example sample sheet {sheet.name}: {exc}")
    module_text="\n".join(path.read_text(encoding="utf-8") for path in (ROOT/"modules/local").rglob("*.nf"))
    for container in re.findall(r"container\s+['\"]([^'\"]+)", module_text):
        if container.endswith(":latest") or ":latest@" in container: errors.append(f"Floating container tag: {container}")
        if ":" not in container: errors.append(f"Unversioned container reference: {container}")
    ignored=(ROOT/".gitignore").read_text(encoding="utf-8")
    for pattern in (".nextflow*", "work/", "results/", ".integration-fixtures/", ".block*_results/", "validation/data/", "validation/work/", "validation/results/", "__pycache__/", "*.pyc", "*.bam", "*.cram", "*.fastq", "*.fastq.gz"):
        if pattern not in ignored: errors.append(f".gitignore missing required protection: {pattern}")
    for path in files():
        relative=path.relative_to(ROOT)
        if path.stat().st_size > 5*1024*1024: errors.append(f"Unexpected source-tree file over 5 MiB: {relative}")
        if path.suffix not in TEXT_SUFFIXES and path.name not in {"Dockerfile", ".gitignore"}: continue
        try: text=path.read_text(encoding="utf-8")
        except UnicodeDecodeError: continue
        if relative.as_posix() != "bin/release_readiness.py":
            if re.search(r"[A-Za-z]:\\Users\\[^\\]+|/home/[^/\s]+|/mnt/[a-z]/Users/", text): errors.append(f"Developer-machine path found in {relative}")
            for pattern in PROHIBITED:
                if re.search(pattern, text, re.IGNORECASE): errors.append(f"Prohibited repository-facing phrase in {relative}: {pattern}")
        if path.suffix == ".md":
            for resolved,target in local_links(path,text):
                if not resolved.exists(): errors.append(f"Broken local documentation link in {relative}: {target}")
        if path.suffix in {".yaml", ".yml", ".cff"}:
            if "\t" in text: errors.append(f"Tab-indented YAML-like file: {relative}")
            if path.suffix != ".cff" and not any(":" in line for line in text.splitlines() if line.strip() and not line.lstrip().startswith("#")):
                errors.append(f"Malformed YAML-like file: {relative}")
    if errors:
        print("Release-readiness checks failed:")
        for item in sorted(set(errors)): print(f"- {item}")
        return 1
    print("Release-readiness checks passed.")
    for item in warnings: print(f"WARNING: {item}")
    return 0


if __name__ == "__main__": sys.exit(main())
