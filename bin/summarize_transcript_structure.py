#!/usr/bin/env python3
"""Validate transcript structure and summarize SQANTI3 classifications."""

import argparse
import csv
import gzip
import json
import statistics
from collections import Counter, defaultdict
from pathlib import Path


def smart_open(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def parse_attributes(text: str) -> dict[str, str]:
    attributes: dict[str, str] = {}
    for part in text.strip().strip(";").split(";"):
        part = part.strip()
        if not part:
            continue
        pieces = part.split(None, 1)
        if len(pieces) != 2:
            raise ValueError(f"malformed GTF attribute: {part!r}")
        attributes[pieces[0]] = pieces[1].strip().strip('"')
    return attributes


def parse_gtf(path: Path) -> dict[str, dict]:
    transcripts: dict[str, dict] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {line_number} does not have 9 columns")
            chrom, _source, feature, start, end, _score, strand, _frame, attr_text = fields
            if strand not in {"+", "-"}:
                raise ValueError(f"GTF line {line_number} has invalid strand {strand!r}")
            try:
                start_i, end_i = int(start), int(end)
            except ValueError as exc:
                raise ValueError(f"GTF line {line_number} has non-integer coordinates") from exc
            if start_i < 1 or end_i < start_i:
                raise ValueError(f"GTF line {line_number} has invalid coordinates")
            attrs = parse_attributes(attr_text)
            transcript_id = attrs.get("transcript_id")
            if feature in {"transcript", "exon"} and not transcript_id:
                raise ValueError(f"GTF line {line_number} lacks transcript_id")
            if not transcript_id:
                continue
            record = transcripts.setdefault(
                transcript_id,
                {
                    "transcript_id": transcript_id,
                    "gene_id": attrs.get("gene_id", "NA"),
                    "chrom": chrom,
                    "strand": strand,
                },
            )
            if record["chrom"] != chrom or record["strand"] != strand:
                raise ValueError(f"transcript {transcript_id!r} spans inconsistent loci")
            if feature == "exon":
                exons[transcript_id].append((start_i, end_i))

    for transcript_id in transcripts:
        if not exons[transcript_id]:
            raise ValueError(f"transcript {transcript_id!r} has no exon records")
    return {tid: {**record, "exons": sorted(exons[tid])} for tid, record in transcripts.items()}


def load_classifications(path: Path) -> list[dict[str, str]]:
    if path.stat().st_size == 0:
        return []
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        required = {"isoform", "chrom", "strand", "length", "exons", "structural_category"}
        missing = required.difference(reader.fieldnames or [])
        if missing:
            raise ValueError(f"SQANTI3 classification lacks columns: {', '.join(sorted(missing))}")
        return list(reader)


def load_support(path: Path) -> dict[str, float]:
    support: dict[str, float] = {}
    with smart_open(path) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if not rows:
        return support
    for row in rows[1:]:
        if len(row) < 2:
            continue
        try:
            support[row[0]] = float(row[1])
        except ValueError:
            continue
    return support


def numeric_summary(values: list[float]) -> dict[str, float | int | None]:
    if not values:
        return {"minimum": None, "median": None, "mean": None, "maximum": None}
    return {
        "minimum": min(values),
        "median": statistics.median(values),
        "mean": statistics.fmean(values),
        "maximum": max(values),
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--support-counts", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--structure-output", required=True, type=Path)
    parser.add_argument("--json-output", required=True, type=Path)
    parser.add_argument("--markdown-output", required=True, type=Path)
    args = parser.parse_args()

    transcripts = parse_gtf(args.gtf)
    classifications = load_classifications(args.classification)
    support = load_support(args.support_counts)

    classification_by_id = {row["isoform"]: row for row in classifications}
    category_counts = Counter(row["structural_category"] for row in classifications)
    exon_counts = [len(record["exons"]) for record in transcripts.values()]
    lengths = [sum(end - start + 1 for start, end in record["exons"]) for record in transcripts.values()]
    support_values = [support[tid] for tid in transcripts if tid in support]

    with args.structure_output.open("w", newline="") as handle:
        fieldnames = [
            "transcript_id", "gene_id", "chrom", "strand", "transcript_start",
            "transcript_end", "transcript_length", "exon_count", "exon_class",
            "exon_coordinates", "intron_chain", "structural_category",
            "associated_gene", "associated_transcript", "supporting_read_count",
        ]
        writer = csv.DictWriter(handle, fieldnames=fieldnames, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for transcript_id in sorted(transcripts):
            record = transcripts[transcript_id]
            exon_list = record["exons"]
            introns = [(exon_list[i][1] + 1, exon_list[i + 1][0] - 1) for i in range(len(exon_list) - 1)]
            classification = classification_by_id.get(transcript_id, {})
            writer.writerow({
                "transcript_id": transcript_id,
                "gene_id": record["gene_id"],
                "chrom": record["chrom"],
                "strand": record["strand"],
                "transcript_start": min(start for start, _end in exon_list),
                "transcript_end": max(end for _start, end in exon_list),
                "transcript_length": sum(end - start + 1 for start, end in exon_list),
                "exon_count": len(exon_list),
                "exon_class": "single-exon" if len(exon_list) == 1 else "multi-exon",
                "exon_coordinates": ",".join(f"{start}-{end}" for start, end in exon_list),
                "intron_chain": ",".join(f"{start}-{end}" for start, end in introns) or "NA",
                "structural_category": classification.get("structural_category", "unclassified"),
                "associated_gene": classification.get("associated_gene", "NA"),
                "associated_transcript": classification.get("associated_transcript", "NA"),
                "supporting_read_count": support.get(transcript_id, "NA"),
            })

    reference_match = category_counts["full-splice_match"] + category_counts["incomplete-splice_match"]
    novel_splice = category_counts["novel_in_catalog"] + category_counts["novel_not_in_catalog"]
    summary = {
        "sample": args.sample,
        "total_transcript_models": len(transcripts),
        "classified_transcript_models": len(classifications),
        "sqanti3_structural_category_counts": dict(sorted(category_counts.items())),
        "reference_matching_fsm_or_ism": reference_match,
        "novel_splice_structure_nic_or_nnc": novel_splice,
        "single_exon_transcripts": sum(count == 1 for count in exon_counts),
        "multi_exon_transcripts": sum(count > 1 for count in exon_counts),
        "transcript_length": numeric_summary([float(value) for value in lengths]),
        "exon_count": numeric_summary([float(value) for value in exon_counts]),
        "supporting_read_count": numeric_summary(support_values),
    }
    args.json_output.write_text(json.dumps(summary, indent=2) + "\n")

    lines = [
        f"# Transcript structural summary: {args.sample}",
        "",
        f"- Total transcript models: {len(transcripts)}",
        f"- SQANTI3-classified models: {len(classifications)}",
        f"- Single-exon models: {summary['single_exon_transcripts']}",
        f"- Multi-exon models: {summary['multi_exon_transcripts']}",
        "",
        "## SQANTI3 structural categories",
        "",
    ]
    if category_counts:
        lines.extend(f"- `{category}`: {count}" for category, count in sorted(category_counts.items()))
    else:
        lines.append("No transcript models were available for classification.")
    args.markdown_output.write_text("\n".join(lines) + "\n")


if __name__ == "__main__":
    main()
