#!/usr/bin/env python3
"""Convert a caller GTF and common SQANTI3 output to the Block 3 structure schema."""

import argparse
import csv
import gzip
import json
from collections import Counter, defaultdict
from pathlib import Path


def open_text(path: Path):
    return gzip.open(path, "rt") if path.suffix == ".gz" else path.open()


def attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.strip().strip(";").split(";"):
        item = item.strip()
        if not item:
            continue
        key, sep, value = item.partition(" ")
        if not sep:
            raise ValueError(f"malformed GTF attribute {item!r}")
        result[key] = value.strip().strip('"')
    return result


def parse_gtf(path: Path) -> dict[str, dict]:
    records: dict[str, dict] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with path.open() as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"GTF line {number} does not have 9 columns")
            chrom, _source, feature, start, end, _score, strand, _frame, raw = fields
            if feature not in {"transcript", "exon"}:
                continue
            attrs = attributes(raw)
            tid = attrs.get("transcript_id")
            if not tid:
                raise ValueError(f"GTF line {number} lacks transcript_id")
            if strand not in {"+", "-"}:
                raise ValueError(f"GTF line {number} has invalid strand")
            start_i, end_i = int(start), int(end)
            if start_i < 1 or end_i < start_i:
                raise ValueError(f"GTF line {number} has invalid coordinates")
            current = records.setdefault(tid, {
                "transcript_id": tid, "gene_id": attrs.get("gene_id", "NA"),
                "chrom": chrom, "strand": strand,
            })
            if (current["chrom"], current["strand"]) != (chrom, strand):
                raise ValueError(f"transcript {tid!r} spans inconsistent loci")
            if feature == "exon":
                exons[tid].append((start_i, end_i))
    for tid, record in records.items():
        if not exons[tid]:
            raise ValueError(f"transcript {tid!r} has no exon records")
        ordered = sorted(set(exons[tid]))
        for left, right in zip(ordered, ordered[1:]):
            if left[1] >= right[0]:
                raise ValueError(f"transcript {tid!r} has overlapping exons")
        record["exons"] = ordered
    return records


def classifications(path: Path) -> dict[str, dict]:
    if not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "isoform" not in reader.fieldnames:
            raise ValueError("SQANTI3 classification lacks isoform column")
        return {row["isoform"]: row for row in reader}


def sqanti_junction_metadata(path: Path | None) -> dict[str, list[dict]]:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return {}
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        if not reader.fieldnames or "isoform" not in reader.fieldnames:
            return {}
        result: dict[str, list[dict]] = defaultdict(list)
        for row in reader:
            result[row["isoform"]].append(row)
    for tid in result:
        result[tid].sort(key=lambda row: int(row.get("junction_number", "0") or 0))
    return result


def support_counts(path: Path | None, kind: str) -> dict[str, str]:
    if path is None or not path.exists() or path.stat().st_size == 0:
        return {}
    counts: Counter[str] = Counter()
    values: dict[str, str] = {}
    with open_text(path) as handle:
        rows = [line.rstrip("\n").split("\t") for line in handle if line.strip()]
    if len(rows) < 2:
        return {}
    header = [item.lower().lstrip("#") for item in rows[0]]
    id_names = ("transcript_id", "isoform", "feature_id", "id")
    count_names = ("count", "read_count", "support", "num_reads")
    id_index = next((header.index(name) for name in id_names if name in header), None)
    count_index = next((header.index(name) for name in count_names if name in header), None)
    if kind in {"read_map", "auto"} and id_index is None and len(header) >= 2:
        id_index = 1
    for row in rows[1:]:
        if id_index is None or id_index >= len(row):
            continue
        tid = row[id_index]
        if not tid:
            continue
        if count_index is not None and count_index < len(row):
            try:
                values[tid] = str(float(row[count_index]))
            except ValueError:
                pass
        else:
            counts[tid] += 1
    values.update({tid: str(value) for tid, value in counts.items()})
    return values


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", required=True, choices=("ont", "pacbio"))
    parser.add_argument("--condition", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--caller", required=True, choices=("isoquant", "flair", "bambu"))
    parser.add_argument("--gtf", required=True, type=Path)
    parser.add_argument("--classification", required=True, type=Path)
    parser.add_argument("--sqanti-junctions", type=Path)
    parser.add_argument("--support", type=Path)
    parser.add_argument("--support-kind", default="auto", choices=("auto", "read_map", "count"))
    parser.add_argument("--transcripts", required=True, type=Path)
    parser.add_argument("--junctions", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()

    models = parse_gtf(args.gtf)
    sqanti = classifications(args.classification)
    sqanti_junctions = sqanti_junction_metadata(args.sqanti_junctions)
    support = support_counts(args.support, args.support_kind)
    fields = [
        "sample", "platform", "condition", "replicate", "caller", "caller_transcript_id",
        "chrom", "strand", "transcript_start", "transcript_end", "tss", "tes",
        "transcript_length", "exon_count", "exon_class", "exon_coordinates", "intron_chain",
        "gene_assignment", "caller_gene_id", "sqanti_category", "sqanti_associated_transcript",
        "support_value", "support_semantics",
    ]
    junction_rows = []
    categories = Counter()
    chains = set()
    junction_set = set()
    with args.transcripts.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        for tid in sorted(models):
            record = models[tid]
            exon_list = record["exons"]
            introns = [(exon_list[i][1] + 1, exon_list[i + 1][0] - 1) for i in range(len(exon_list) - 1)]
            chain = ",".join(f"{start}-{end}" for start, end in introns) or "NA"
            start, end = exon_list[0][0], exon_list[-1][1]
            category = sqanti.get(tid, {}).get("structural_category", "unclassified") or "unclassified"
            associated_gene = sqanti.get(tid, {}).get("associated_gene", "NA") or "NA"
            categories[category] += 1
            if introns:
                chains.add((record["chrom"], record["strand"], chain))
            per_transcript_junctions = sqanti_junctions.get(tid, [])
            for intron_index, (donor, acceptor) in enumerate(introns):
                sqanti_junction = per_transcript_junctions[intron_index] if len(per_transcript_junctions) == len(introns) else {}
                key = (record["chrom"], record["strand"], donor, acceptor)
                junction_set.add(key)
                junction_rows.append({
                    "sample": args.sample, "platform": args.platform, "caller": args.caller,
                    "caller_transcript_id": tid, "chrom": record["chrom"], "strand": record["strand"],
                    "intron_start": donor, "intron_end": acceptor,
                    "canonical_status": sqanti_junction.get("canonical", "NA") or "NA",
                    "reference_status": sqanti_junction.get("junction_category", "NA") or "NA",
                    "support_value": support.get(tid, "NA"),
                    "support_semantics": {
                        "isoquant": "IsoQuant transcript-model read assignment/count",
                        "flair": "FLAIR read-to-isoform assignment count",
                        "bambu": "Bambu caller-native transcript abundance/support",
                    }[args.caller],
                })
            writer.writerow({
                "sample": args.sample, "platform": args.platform, "condition": args.condition,
                "replicate": args.replicate, "caller": args.caller, "caller_transcript_id": tid,
                "chrom": record["chrom"], "strand": record["strand"],
                "transcript_start": start, "transcript_end": end,
                "tss": start if record["strand"] == "+" else end,
                "tes": end if record["strand"] == "+" else start,
                "transcript_length": sum(b - a + 1 for a, b in exon_list),
                "exon_count": len(exon_list),
                "exon_class": "single-exon" if len(exon_list) == 1 else "multi-exon",
                "exon_coordinates": ",".join(f"{a}-{b}" for a, b in exon_list),
                "intron_chain": chain,
                "gene_assignment": associated_gene if associated_gene != "NA" else record["gene_id"],
                "caller_gene_id": record["gene_id"], "sqanti_category": category,
                "sqanti_associated_transcript": sqanti.get(tid, {}).get("associated_transcript", "NA") or "NA",
                "support_value": support.get(tid, "NA"),
                "support_semantics": {
                    "isoquant": "IsoQuant transcript-model read assignment/count",
                    "flair": "FLAIR read-to-isoform assignment count",
                    "bambu": "Bambu caller-native transcript abundance/support",
                }[args.caller],
            })

    junction_fields = ["sample", "platform", "caller", "caller_transcript_id", "chrom", "strand",
                       "intron_start", "intron_end", "canonical_status", "reference_status", "support_value",
                       "support_semantics"]
    with args.junctions.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=junction_fields, delimiter="\t", lineterminator="\n")
        writer.writeheader()
        writer.writerows(sorted(junction_rows, key=lambda row: (row["chrom"], row["strand"], row["intron_start"], row["intron_end"], row["caller_transcript_id"])))

    result = {
        "sample": args.sample, "platform": args.platform, "condition": args.condition,
        "replicate": args.replicate, "caller": args.caller,
        "transcripts": len(models),
        "multi_exon_transcripts": sum(len(item["exons"]) > 1 for item in models.values()),
        "single_exon_transcripts": sum(len(item["exons"]) == 1 for item in models.values()),
        "unique_intron_chains": len(chains), "splice_junctions": len(junction_set),
        "sqanti3_category_counts": dict(sorted(categories.items())),
        "support_semantics": {
            "isoquant": "caller-native IsoQuant read assignment/count; not equated to other callers",
            "flair": "number of rows in the FLAIR read-to-isoform map; not equated to other callers",
            "bambu": "caller-native Bambu abundance/support; not equated to other callers",
        }[args.caller],
    }
    args.summary.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
