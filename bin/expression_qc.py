#!/usr/bin/env python3
"""Descriptive abundance QC; never excludes samples automatically."""

import argparse
import csv
import json
import math
from pathlib import Path


def matrix(path: Path):
    with path.open(newline="") as handle:
        reader = csv.reader(handle, delimiter="\t")
        header = next(reader); samples = header[1:]
        rows = [(row[0], [float(x) for x in row[1:]]) for row in reader]
    return samples, rows


def pearson(left, right):
    if not left or len(left) != len(right): return None
    ml, mr = sum(left) / len(left), sum(right) / len(right)
    numerator = sum((a - ml) * (b - mr) for a, b in zip(left, right))
    dl = sum((a - ml) ** 2 for a in left); dr = sum((b - mr) ** 2 for b in right)
    return None if dl == 0 or dr == 0 else numerator / math.sqrt(dl * dr)


def write_bar_svg(path: Path, labels, values, title):
    width, height, margin = 700, 400, 70
    maximum = max(values, default=0) or 1
    bar_width = max(10, (width - 2 * margin) / max(len(values), 1) * 0.65)
    pieces = [f'<svg xmlns="http://www.w3.org/2000/svg" width="{width}" height="{height}">',
              f'<text x="{width/2}" y="24" text-anchor="middle">{title}</text>']
    for index, (label, value) in enumerate(zip(labels, values)):
        x = margin + index * (width - 2 * margin) / max(len(values), 1)
        h = value / maximum * (height - 2 * margin)
        pieces += [f'<rect x="{x}" y="{height-margin-h}" width="{bar_width}" height="{h}" fill="#3977a8"/>',
                   f'<text x="{x+bar_width/2}" y="{height-margin+16}" text-anchor="middle" font-size="10">{label}</text>',
                   f'<text x="{x+bar_width/2}" y="{height-margin-h-4}" text-anchor="middle" font-size="10">{value:g}</text>']
    pieces.append('</svg>')
    path.write_text("\n".join(pieces) + "\n")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--dataset", required=True, type=Path)
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args(); args.output.mkdir(parents=True, exist_ok=True)
    samples, genes = matrix(args.dataset / "gene_counts.tsv")
    _, transcripts = matrix(args.dataset / "transcript_counts.tsv")
    gene_vectors = [[math.log1p(row[1][i]) for row in genes] for i in range(len(samples))]
    transcript_vectors = [[math.log1p(row[1][i]) for row in transcripts] for i in range(len(samples))]
    rows = []
    for i, sample in enumerate(samples):
        gene_counts = [row[1][i] for row in genes]; tx_counts = [row[1][i] for row in transcripts]
        rows.append({"sample": sample, "library_count_total": int(round(sum(gene_counts))),
                     "transcript_count_total": int(round(sum(tx_counts))),
                     "detected_genes": sum(x > 0 for x in gene_counts),
                     "detected_transcripts": sum(x > 0 for x in tx_counts)})
    with (args.output / "sample_abundance_qc.tsv").open("w", newline="") as handle:
        fields = ["sample", "library_count_total", "transcript_count_total", "detected_genes", "detected_transcripts"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n"); writer.writeheader(); writer.writerows(rows)
    with (args.output / "sample_correlation.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(["sample", *samples])
        for i, sample in enumerate(samples):
            writer.writerow([sample, *[("NA" if (value := pearson(transcript_vectors[i], transcript_vectors[j])) is None else f"{value:.8g}") for j in range(len(samples))]])
    pca_status = "NOT_RUN_TOO_FEW_SAMPLES" if len(samples) < 3 else "NOT_RUN_NO_VARIABLE_FEATURES"
    try:
        import numpy as np
        array = np.array(transcript_vectors, dtype=float)
        if len(samples) >= 3 and array.shape[1] >= 2 and np.any(np.var(array, axis=0) > 0):
            centered = array - array.mean(axis=0)
            u, s, _ = np.linalg.svd(centered, full_matrices=False)
            with (args.output / "pca_scores.tsv").open("w", newline="") as handle:
                writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(["sample", "PC1", "PC2"])
                for index, sample in enumerate(samples): writer.writerow([sample, u[index, 0] * s[0], u[index, 1] * s[1] if len(s) > 1 else 0])
            pca_status = "COMPLETE"
    except ImportError:
        pca_status = "NOT_RUN_NUMPY_UNAVAILABLE"
    (args.output / "expression_qc_status.json").write_text(json.dumps({"status": "COMPLETE", "pca_status": pca_status,
        "sample_exclusion_performed": False, "interpretation": "descriptive_qc_only"}, indent=2) + "\n")
    write_bar_svg(args.output / "library_totals.svg", samples, [x["library_count_total"] for x in rows], "Raw assigned gene-count totals")
    write_bar_svg(args.output / "detected_transcripts.svg", samples, [x["detected_transcripts"] for x in rows], "Detected production transcripts")


if __name__ == "__main__": main()
