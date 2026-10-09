#!/usr/bin/env python3
"""Compare standardized transcript structures without treating IDs as identity."""

import argparse
import csv
import json
from collections import Counter, defaultdict
from itertools import combinations
from pathlib import Path

CALLERS = ("isoquant", "flair", "bambu")
PAIR_NAMES = (("isoquant", "flair"), ("isoquant", "bambu"), ("flair", "bambu"))


def load(path: Path) -> list[dict]:
    with path.open(newline="") as handle:
        rows = list(csv.DictReader(handle, delimiter="\t"))
    for row in rows:
        for field in ("transcript_start", "transcript_end", "tss", "tes", "exon_count"):
            row[field] = int(row[field])
    return rows


def exact_key(row: dict) -> tuple:
    return (row["chrom"], row["strand"], row["intron_chain"], row["transcript_start"], row["transcript_end"])


def chain_key(row: dict) -> tuple:
    return (row["chrom"], row["strand"], row["intron_chain"])


def junctions(row: dict) -> set[tuple]:
    if row["intron_chain"] == "NA":
        return set()
    return {(row["chrom"], row["strand"], *map(int, value.split("-"))) for value in row["intron_chain"].split(",")}


def jaccard(left: set, right: set) -> dict:
    intersection = len(left & right)
    union = len(left | right)
    if not union:
        return {"intersection": 0, "union": 0, "jaccard": None, "status": "both_empty"}
    return {"intersection": intersection, "union": union, "jaccard": intersection / union,
            "status": "one_empty" if not left or not right else "defined"}


def overlap(a: dict, b: dict) -> int:
    if (a["chrom"], a["strand"]) != (b["chrom"], b["strand"]):
        return 0
    return max(0, min(a["transcript_end"], b["transcript_end"]) - max(a["transcript_start"], b["transcript_start"]) + 1)


def end_bin(value: int, bins: list[int]) -> str:
    for limit in bins:
        if value <= limit:
            return "exact" if limit == 0 else f"<={limit}bp"
    return f">{bins[-1]}bp"


def category_group(category: str) -> str:
    mapping = {
        "full-splice_match": "FSM", "incomplete-splice_match": "ISM",
        "novel_in_catalog": "NIC", "novel_not_in_catalog": "NNC",
        "antisense": "antisense", "intergenic": "intergenic",
    }
    return mapping.get(category, category or "unclassified")


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def main() -> None:
    parser = argparse.ArgumentParser()
    for caller in CALLERS:
        parser.add_argument(f"--{caller}", required=True, type=Path)
    parser.add_argument("--sample", required=True)
    parser.add_argument("--platform", required=True)
    parser.add_argument("--condition", required=True)
    parser.add_argument("--replicate", required=True)
    parser.add_argument("--end-bins", default="0,10,25,50")
    parser.add_argument("--single-exon-overlap", default=0.5, type=float)
    parser.add_argument("--membership", required=True, type=Path)
    parser.add_argument("--pairwise", required=True, type=Path)
    parser.add_argument("--ends", required=True, type=Path)
    parser.add_argument("--partial-junctions", required=True, type=Path)
    parser.add_argument("--gene-overlap", required=True, type=Path)
    parser.add_argument("--summary", required=True, type=Path)
    args = parser.parse_args()
    bins = sorted(set(int(value) for value in args.end_bins.split(",")))
    if not bins or bins[0] != 0 or any(value < 0 for value in bins):
        raise ValueError("end bins must be non-negative and include 0")
    if not 0 < args.single_exon_overlap <= 1:
        raise ValueError("single-exon overlap must be in (0, 1]")

    data = {caller: load(getattr(args, caller)) for caller in CALLERS}
    multi = {caller: [row for row in rows if row["exon_count"] > 1] for caller, rows in data.items()}
    single = {caller: [row for row in rows if row["exon_count"] == 1] for caller, rows in data.items()}
    exact_sets = {caller: {exact_key(row) for row in multi[caller]} for caller in CALLERS}
    chain_sets = {caller: {chain_key(row) for row in multi[caller]} for caller in CALLERS}
    junction_sets = {caller: set().union(*(junctions(row) for row in multi[caller])) if multi[caller] else set() for caller in CALLERS}

    # Membership uses exact multi-exon structure. Single-exon models form groups by reciprocal genomic overlap.
    groups = []
    all_exact = sorted(set().union(*exact_sets.values()))
    for key in all_exact:
        members = {caller: sorted(row["caller_transcript_id"] for row in multi[caller] if exact_key(row) == key) for caller in CALLERS}
        exemplar = next(row for caller in CALLERS for row in multi[caller] if exact_key(row) == key)
        groups.append(("multi-exon exact", exemplar, members))
    single_nodes = [(caller, row) for caller in CALLERS for row in single[caller]]
    parent = list(range(len(single_nodes)))
    def find(x):
        while parent[x] != x:
            parent[x] = parent[parent[x]]
            x = parent[x]
        return x
    def union(a, b):
        a, b = find(a), find(b)
        if a != b:
            parent[b] = a
    for (i, (ca, a)), (j, (cb, b)) in combinations(enumerate(single_nodes), 2):
        if ca == cb:
            continue
        shared = overlap(a, b)
        if shared and shared / (a["transcript_end"] - a["transcript_start"] + 1) >= args.single_exon_overlap and shared / (b["transcript_end"] - b["transcript_start"] + 1) >= args.single_exon_overlap:
            union(i, j)
    components = defaultdict(list)
    for index, node in enumerate(single_nodes):
        components[find(index)].append(node)
    for nodes in sorted(components.values(), key=lambda n: (n[0][1]["chrom"], n[0][1]["transcript_start"], n[0][0])):
        exemplar = nodes[0][1]
        members = {caller: sorted(row["caller_transcript_id"] for source, row in nodes if source == caller) for caller in CALLERS}
        groups.append(("single-exon reciprocal overlap", exemplar, members))

    membership_rows = []
    membership_counts = Counter()
    category_membership = defaultdict(Counter)
    for number, (basis, exemplar, members) in enumerate(groups, 1):
        present = tuple(caller for caller in CALLERS if members[caller])
        label = "+".join(present) if len(present) > 1 else f"{present[0]}_only"
        membership_counts[label] += 1
        categories = sorted({category_group(row["sqanti_category"]) for caller in CALLERS for row in data[caller] if row["caller_transcript_id"] in members[caller]})
        for category in categories or ["unclassified"]:
            category_membership[category][label] += 1
        membership_rows.append({
            "sample": args.sample, "platform": args.platform, "structure_group": f"SG{number:06d}",
            "matching_basis": basis, "chrom": exemplar["chrom"], "strand": exemplar["strand"],
            "isoquant_transcript_ids": ",".join(members["isoquant"]) or "NA",
            "flair_transcript_ids": ",".join(members["flair"]) or "NA",
            "bambu_transcript_ids": ",".join(members["bambu"]) or "NA",
            "membership": label, "sqanti_category_groups": ",".join(categories) or "unclassified",
        })

    pairwise_rows, end_rows, partial_rows, gene_rows = [], [], [], []
    pairwise_json = {}
    for left, right in PAIR_NAMES:
        pair = f"{left}_vs_{right}"
        metrics = {
            "exact_multi_exon_structures": jaccard(exact_sets[left], exact_sets[right]),
            "intron_chains": jaccard(chain_sets[left], chain_sets[right]),
            "splice_junctions": jaccard(junction_sets[left], junction_sets[right]),
        }
        # Same exact coordinates are reported separately for single-exon models; they never enter intron-chain metrics.
        left_single_exact = {(r["chrom"], r["strand"], r["transcript_start"], r["transcript_end"]) for r in single[left]}
        right_single_exact = {(r["chrom"], r["strand"], r["transcript_start"], r["transcript_end"]) for r in single[right]}
        metrics["exact_single_exon_coordinates"] = jaccard(left_single_exact, right_single_exact)
        pairwise_json[pair] = metrics
        for level, result in metrics.items():
            pairwise_rows.append({"sample": args.sample, "platform": args.platform, "caller_pair": pair, "level": level, **result})
        for a in multi[left]:
            for b in multi[right]:
                if chain_key(a) == chain_key(b):
                    dtss, dtes = b["tss"] - a["tss"], b["tes"] - a["tes"]
                    end_rows.append({"sample": args.sample, "platform": args.platform, "caller_pair": pair,
                                     "left_transcript_id": a["caller_transcript_id"], "right_transcript_id": b["caller_transcript_id"],
                                     "chrom": a["chrom"], "strand": a["strand"], "tss_difference": dtss,
                                     "tes_difference": dtes, "absolute_tss_difference": abs(dtss), "absolute_tes_difference": abs(dtes),
                                     "tss_bin": end_bin(abs(dtss), bins), "tes_bin": end_bin(abs(dtes), bins)})
                else:
                    shared = junctions(a) & junctions(b)
                    if shared:
                        partial_rows.append({"sample": args.sample, "platform": args.platform, "caller_pair": pair,
                                             "left_transcript_id": a["caller_transcript_id"], "right_transcript_id": b["caller_transcript_id"],
                                             "shared_junctions": len(shared), "left_junctions": len(junctions(a)), "right_junctions": len(junctions(b))})
        left_gene = defaultdict(list); right_gene = defaultdict(list)
        for row in data[left]:
            if row["gene_assignment"] != "NA": left_gene[row["gene_assignment"]].append(row["caller_transcript_id"])
        for row in data[right]:
            if row["gene_assignment"] != "NA": right_gene[row["gene_assignment"]].append(row["caller_transcript_id"])
        for gene in sorted(left_gene.keys() & right_gene.keys()):
            gene_rows.append({"sample": args.sample, "platform": args.platform, "caller_pair": pair, "gene_assignment": gene,
                              "left_transcript_ids": ",".join(sorted(left_gene[gene])), "right_transcript_ids": ",".join(sorted(right_gene[gene]))})

    write_tsv(args.membership, ["sample", "platform", "structure_group", "matching_basis", "chrom", "strand",
                                      "isoquant_transcript_ids", "flair_transcript_ids", "bambu_transcript_ids", "membership", "sqanti_category_groups"], membership_rows)
    write_tsv(args.pairwise, ["sample", "platform", "caller_pair", "level", "intersection", "union", "jaccard", "status"], pairwise_rows)
    write_tsv(args.ends, ["sample", "platform", "caller_pair", "left_transcript_id", "right_transcript_id", "chrom", "strand",
                                "tss_difference", "tes_difference", "absolute_tss_difference", "absolute_tes_difference", "tss_bin", "tes_bin"], end_rows)
    write_tsv(args.partial_junctions, ["sample", "platform", "caller_pair", "left_transcript_id", "right_transcript_id",
                                             "shared_junctions", "left_junctions", "right_junctions"], partial_rows)
    write_tsv(args.gene_overlap, ["sample", "platform", "caller_pair", "gene_assignment", "left_transcript_ids", "right_transcript_ids"], gene_rows)

    summary = {
        "contract_version": "3.0", "sample": args.sample, "platform": args.platform,
        "condition": args.condition, "replicate": args.replicate,
        "identity_policy": {"multi_exon": "chromosome + strand + ordered intron chain + exact transcript ends",
                            "single_exon": f"chromosome + strand + reciprocal genomic overlap >= {args.single_exon_overlap}",
                            "caller_transcript_ids_used_for_identity": False},
        "end_difference_bins_bp": bins,
        "caller_counts": {caller: {"transcripts": len(data[caller]), "multi_exon": len(multi[caller]), "single_exon": len(single[caller]),
                                    "unique_intron_chains": len(chain_sets[caller]), "splice_junctions": len(junction_sets[caller]),
                                    "sqanti3_category_counts": dict(sorted(Counter(row["sqanti_category"] for row in data[caller]).items()))}
                          for caller in CALLERS},
        "membership_counts": dict(sorted(membership_counts.items())),
        "membership_by_sqanti_category": {category: dict(sorted(counts.items())) for category, counts in sorted(category_membership.items())},
        "pairwise": pairwise_json,
        "intron_chain_matched_end_pairs": len(end_rows), "partial_junction_pairs": len(partial_rows),
        "consensus_transcriptome_created": False,
    }
    args.summary.write_text(json.dumps(summary, indent=2) + "\n")


if __name__ == "__main__":
    main()
