#!/usr/bin/env python3
"""Summarize observed transcript boundaries and configurable end clusters."""

import argparse
import csv
import json
from collections import defaultdict
from pathlib import Path


def read_matrix(path):
    with path.open(newline="") as handle:
        reader = csv.DictReader(handle, delimiter="\t")
        return reader.fieldnames[1:], {row["feature_id"]: row for row in reader}


def cluster(values, window):
    result, cluster_id, anchor = {}, 0, None
    for value in sorted(set(values)):
        if anchor is None or value - anchor > window:
            cluster_id += 1; anchor = value
        result[value] = cluster_id
    return result


def main():
    p = argparse.ArgumentParser(); p.add_argument("--dataset", type=Path, required=True); p.add_argument("--window", type=int, default=50); p.add_argument("--output", type=Path, required=True)
    a = p.parse_args(); a.output.mkdir(parents=True, exist_ok=True)
    if a.window < 0: raise ValueError("Boundary clustering window must be non-negative")
    samples, tpm = read_matrix(a.dataset / "transcript_tpm.tsv")
    with (a.dataset / "transcript_metadata.tsv").open(newline="") as handle: models = list(csv.DictReader(handle, delimiter="\t"))
    by_gene = defaultdict(list)
    for row in models:
        start, end = int(row["transcript_start"]), int(row["transcript_end"])
        row["observed_transcript_start"] = start; row["observed_transcript_end"] = end
        row["putative_tss"] = start if row["strand"] == "+" else end
        row["putative_tes"] = end if row["strand"] == "+" else start
        row["sample_support"] = sum(float(tpm.get(row["transcript_id"], {}).get(s, 0)) > 0 for s in samples)
        row["abundance_sum_tpm"] = sum(float(tpm.get(row["transcript_id"], {}).get(s, 0)) for s in samples)
        by_gene[row.get("gene_id") or row.get("associated_gene") or "NA"].append(row)
    for gene, rows in by_gene.items():
        tss_clusters = cluster([x["putative_tss"] for x in rows], a.window); tes_clusters = cluster([x["putative_tes"] for x in rows], a.window)
        for row in rows:
            row["tss_cluster"] = f"{gene}:TSS:{tss_clusters[row['putative_tss']]}"
            row["tes_cluster"] = f"{gene}:TES:{tes_clusters[row['putative_tes']]}"
    fields = ["transcript_id", "gene_id", "chrom", "strand", "observed_transcript_start", "observed_transcript_end", "putative_tss", "putative_tes", "transcript_length", "structural_category", "known_novel", "sample_support", "abundance_sum_tpm", "tss_cluster", "tes_cluster"]
    with (a.output / "transcript_boundaries.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n", extrasaction="ignore"); writer.writeheader(); writer.writerows(models)
    with (a.output / "gene_boundary_summary.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n"); writer.writerow(["gene_id", "transcripts", "distinct_observed_starts", "distinct_observed_ends", "tss_clusters", "tes_clusters", "abundance_sum_tpm"])
        for gene, rows in sorted(by_gene.items()): writer.writerow([gene, len(rows), len({x['putative_tss'] for x in rows}), len({x['putative_tes'] for x in rows}), len({x['tss_cluster'] for x in rows}), len({x['tes_cluster'] for x in rows}), sum(x['abundance_sum_tpm'] for x in rows)])
    (a.output / "transcript_ends_status.json").write_text(json.dumps({"status":"COMPLETE", "cluster_window_bp":a.window, "orthogonal_validation":"NOT_PROVIDED", "terminology":"putative TSS/TES"}, indent=2)+"\n")


if __name__ == "__main__": main()
