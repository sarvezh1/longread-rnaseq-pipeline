#!/usr/bin/env python3
"""Deterministic Block 6 truth, reproducibility, depth, and abundance metrics."""

from __future__ import annotations

import argparse
import csv
import gzip
import hashlib
import json
import math
import statistics
from collections import defaultdict
from pathlib import Path


def open_text(path: Path, mode: str = "rt"):
    return gzip.open(path, mode, encoding="utf-8") if path.suffix == ".gz" else path.open(mode, encoding="utf-8", newline="")


def read_tsv(path: Path) -> list[dict[str, str]]:
    with open_text(path) as handle:
        return list(csv.DictReader(handle, delimiter="\t"))


def write_tsv(path: Path, fields: list[str], rows: list[dict]) -> None:
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", encoding="utf-8", newline="") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        writer.writerows(rows)


def file_sha256(path: Path) -> str:
    value=hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda:handle.read(1024*1024),b""): value.update(block)
    return value.hexdigest()


def attributes(text: str) -> dict[str, str]:
    result = {}
    for item in text.strip().strip(";").split(";"):
        key, separator, value = item.strip().partition(" ")
        if separator:
            result[key] = value.strip().strip('"')
    return result


def models_from_gtf(path: Path) -> list[dict]:
    loci: dict[str, tuple[str, str, str]] = {}
    exons: dict[str, list[tuple[int, int]]] = defaultdict(list)
    with open_text(path) as handle:
        for number, line in enumerate(handle, 1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise ValueError(f"{path}: GTF line {number} does not have 9 columns")
            chrom, _source, feature, start, end, _score, strand, _frame, raw = fields
            if feature != "exon":
                continue
            attrs = attributes(raw)
            transcript = attrs.get("transcript_id")
            if not transcript:
                raise ValueError(f"{path}: exon line {number} lacks transcript_id")
            locus = (chrom, strand, attrs.get("gene_id", "NA"))
            if transcript in loci and loci[transcript][:2] != locus[:2]:
                raise ValueError(f"{path}: transcript {transcript} spans inconsistent loci")
            loci[transcript] = locus
            exons[transcript].append((int(start), int(end)))
    result = []
    for transcript in sorted(exons):
        chrom, strand, gene = loci[transcript]
        ordered = sorted(set(exons[transcript]))
        if strand not in {"+", "-"} or any(a < 1 or b < a for a, b in ordered):
            raise ValueError(f"{path}: invalid coordinates or strand for {transcript}")
        if any(left[1] >= right[0] for left, right in zip(ordered, ordered[1:])):
            raise ValueError(f"{path}: overlapping exons for {transcript}")
        introns = tuple((ordered[i][1] + 1, ordered[i + 1][0] - 1) for i in range(len(ordered) - 1))
        start, end = ordered[0][0], ordered[-1][1]
        result.append({"transcript_id": transcript, "gene_id": gene, "chrom": chrom, "strand": strand,
                       "start": start, "end": end, "tss": start if strand == "+" else end,
                       "tes": end if strand == "+" else start, "exons": tuple(ordered), "introns": introns,
                       "exon_count": len(ordered)})
    return result


def models_from_standardized(path: Path) -> list[dict]:
    rows = read_tsv(path)
    result = []
    for row in rows:
        exon_text = row.get("exon_coordinates", "")
        exons = tuple(tuple(map(int, part.split("-"))) for part in exon_text.split(",") if part and part != "NA")
        if not exons:
            start, end = int(row["transcript_start"]), int(row["transcript_end"])
            exons = ((start, end),)
        introns = tuple((exons[i][1] + 1, exons[i + 1][0] - 1) for i in range(len(exons) - 1))
        strand = row["strand"]
        result.append({**row, "transcript_id": row.get("caller_transcript_id", row.get("transcript_id", "NA")),
                       "gene_id": row.get("gene_assignment", row.get("gene_id", "NA")), "chrom": row["chrom"],
                       "strand": strand, "start": exons[0][0], "end": exons[-1][1],
                       "tss": exons[0][0] if strand == "+" else exons[-1][1],
                       "tes": exons[-1][1] if strand == "+" else exons[0][0],
                       "exons": exons, "introns": introns, "exon_count": len(exons)})
    return result


def unique(models: list[dict], key) -> list[dict]:
    observed = {}
    for model in models:
        observed.setdefault(key(model), model)
    return list(observed.values())


def safe_metric(numerator: int, denominator: int):
    return numerator / denominator if denominator else None


def prf(tp: int, fp: int, fn: int) -> dict:
    if tp + fn == 0:
        return {"tp": tp, "fp": fp, "fn": fn, "precision": None, "recall": None, "f1": None, "status": "EMPTY_TRUTH"}
    if tp + fp == 0:
        return {"tp": tp, "fp": fp, "fn": fn, "precision": None, "recall": 0.0, "f1": None, "status": "ZERO_DETECTION"}
    precision, recall = safe_metric(tp, tp + fp), safe_metric(tp, tp + fn)
    f1 = None if precision is None or recall is None or precision + recall == 0 else 2 * precision * recall / (precision + recall)
    return {"tp": tp, "fp": fp, "fn": fn, "precision": precision, "recall": recall, "f1": f1,
            "status": "COMPLETE"}


def bipartite_pairs(predicted: list[dict], truth: list[dict], compatible) -> list[tuple[int, int]]:
    candidates = sorted((score, pi, ti) for pi, p in enumerate(predicted) for ti, t in enumerate(truth)
                        for score in [compatible(p, t)] if score is not None)
    used_p, used_t, pairs = set(), set(), []
    for _score, pi, ti in candidates:
        if pi not in used_p and ti not in used_t:
            used_p.add(pi); used_t.add(ti); pairs.append((pi, ti))
    return pairs


def score_truth(predicted: list[dict], truth: list[dict], tolerance: int) -> tuple[dict, list[dict]]:
    multi_pred = unique([m for m in predicted if m["exon_count"] > 1], lambda m: (m["chrom"], m["strand"], m["introns"], m["start"], m["end"]))
    multi_truth = unique([m for m in truth if m["exon_count"] > 1], lambda m: (m["chrom"], m["strand"], m["introns"], m["start"], m["end"]))
    single_pred = unique([m for m in predicted if m["exon_count"] == 1], lambda m: (m["chrom"], m["strand"], m["start"], m["end"]))
    single_truth = unique([m for m in truth if m["exon_count"] == 1], lambda m: (m["chrom"], m["strand"], m["start"], m["end"]))
    def exact_compatible(p, t):
        if (p["chrom"], p["strand"], p["introns"]) != (t["chrom"], t["strand"], t["introns"]): return None
        if abs(p["tss"] - t["tss"]) > tolerance or abs(p["tes"] - t["tes"]) > tolerance: return None
        return abs(p["tss"] - t["tss"]) + abs(p["tes"] - t["tes"])
    exact_multi_pairs = bipartite_pairs(multi_pred, multi_truth, exact_compatible)
    exact_single_pairs = bipartite_pairs(single_pred, single_truth, exact_compatible)
    chain_pred = unique(multi_pred, lambda m: (m["chrom"], m["strand"], m["introns"]))
    chain_truth = unique(multi_truth, lambda m: (m["chrom"], m["strand"], m["introns"]))
    chain_pairs = bipartite_pairs(chain_pred, chain_truth, lambda p, t: 0 if (p["chrom"], p["strand"], p["introns"]) == (t["chrom"], t["strand"], t["introns"]) else None)
    pred_junctions = {(m["chrom"], m["strand"], *j) for m in multi_pred for j in m["introns"]}
    truth_junctions = {(m["chrom"], m["strand"], *j) for m in multi_truth for j in m["introns"]}
    exact_total = len(exact_multi_pairs) + len(exact_single_pairs)
    metrics = {
        "matching_definition": {"end_tolerance_bp": tolerance, "multi_exon_exact": "chromosome+strand+ordered intron chain+TSS/TES tolerance", "single_exon_exact": "chromosome+strand+TSS/TES tolerance"},
        "exact_all": prf(exact_total, len(multi_pred) + len(single_pred) - exact_total, len(multi_truth) + len(single_truth) - exact_total),
        "exact_multi_exon": prf(len(exact_multi_pairs), len(multi_pred) - len(exact_multi_pairs), len(multi_truth) - len(exact_multi_pairs)),
        "exact_single_exon": prf(len(exact_single_pairs), len(single_pred) - len(exact_single_pairs), len(single_truth) - len(exact_single_pairs)),
        "intron_chain": prf(len(chain_pairs), len(chain_pred) - len(chain_pairs), len(chain_truth) - len(chain_pairs)),
        "splice_junction": prf(len(pred_junctions & truth_junctions), len(pred_junctions - truth_junctions), len(truth_junctions - pred_junctions)),
    }
    matches = []
    for pi, ti in chain_pairs:
        p, t = chain_pred[pi], chain_truth[ti]
        matches.append({"predicted_transcript_id": p["transcript_id"], "truth_transcript_id": t["transcript_id"],
                        "chrom": p["chrom"], "strand": p["strand"], "tss_error_bp": p["tss"] - t["tss"],
                        "tes_error_bp": p["tes"] - t["tes"], "absolute_tss_error_bp": abs(p["tss"] - t["tss"]),
                        "absolute_tes_error_bp": abs(p["tes"] - t["tes"])})
    for key in ("absolute_tss_error_bp", "absolute_tes_error_bp"):
        values = [row[key] for row in matches]
        metrics[key] = {"n": len(values), "median": statistics.median(values) if values else None,
                        "mean": statistics.fmean(values) if values else None,
                        "maximum": max(values) if values else None, "status": "COMPLETE" if values else "NOT_ESTIMABLE_NO_INTRON_CHAIN_MATCH"}
    return metrics, matches


def cmd_truth(args) -> None:
    prediction = models_from_standardized(args.prediction) if args.prediction_format == "standardized" else models_from_gtf(args.prediction)
    truth = models_from_gtf(args.truth)
    if args.truth_contig_prefix:
        prediction = [model for model in prediction if model["chrom"].startswith(args.truth_contig_prefix)]
        truth = [model for model in truth if model["chrom"].startswith(args.truth_contig_prefix)]
    metrics, matches = score_truth(prediction, truth, args.end_tolerance)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps({"schema_version": "6.0", "dataset_id": args.dataset_id, "caller": args.caller,
                                       "platform": args.platform, "truth_type": args.truth_type, **metrics}, indent=2) + "\n", encoding="utf-8")
    write_tsv(args.matches, ["predicted_transcript_id", "truth_transcript_id", "chrom", "strand", "tss_error_bp", "tes_error_bp", "absolute_tss_error_bp", "absolute_tes_error_bp"], matches)


def cmd_truth_inventory(args) -> None:
    models=models_from_gtf(args.truth)
    if args.contig_prefix: models=[model for model in models if model["chrom"].startswith(args.contig_prefix)]
    junctions={(model["chrom"],model["strand"],*junction) for model in models for junction in model["introns"]}
    result={"schema_version":"6.0","truth_id":args.truth_id,"source_file":args.truth.name,"source_bytes":args.truth.stat().st_size,"source_sha256":file_sha256(args.truth),"contig_prefix":args.contig_prefix or None,"transcripts":len(models),"multi_exon_transcripts":sum(model["exon_count"]>1 for model in models),"single_exon_transcripts":sum(model["exon_count"]==1 for model in models),"splice_junctions":len(junctions),"reference_sequences":len({model["chrom"] for model in models})}
    args.output.parent.mkdir(parents=True,exist_ok=True); args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


def fastq_records(path: Path):
    with open_text(path) as handle:
        while True:
            lines = [handle.readline() for _ in range(4)]
            if not lines[0]: return
            if any(line == "" for line in lines) or not lines[0].startswith("@") or not lines[2].startswith("+"):
                raise ValueError(f"Malformed or truncated FASTQ record in {path}")
            yield lines


def cmd_downsample(args) -> None:
    if not 0 < args.fraction <= 1: raise ValueError("fraction must be in (0,1]")
    threshold = int(args.fraction * (1 << 256)); kept = total = 0
    opener = gzip.open if args.output.suffix == ".gz" else open
    args.output.parent.mkdir(parents=True, exist_ok=True)
    with opener(args.output, "wt", encoding="utf-8") as output:
        for record in fastq_records(args.input):
            total += 1
            read_id = record[0][1:].split()[0]
            digest = int(hashlib.sha256(f"{args.seed}\0{read_id}".encode()).hexdigest(), 16)
            if digest < threshold:
                output.writelines(record); kept += 1
    args.manifest.write_text(json.dumps({"schema_version": "6.0", "method": "SHA-256 read-identifier threshold",
                                         "seed": args.seed, "fraction": args.fraction, "input_reads": total,
                                         "retained_reads": kept, "nested_across_fractions": True}, indent=2) + "\n", encoding="utf-8")


def cmd_fastq_stats(args) -> None:
    lengths=[]
    for record in fastq_records(args.input): lengths.append(len(record[1].rstrip("\r\n")))
    ordered=sorted(lengths,reverse=True); halfway=sum(lengths)/2; cumulative=0; n50=None
    for length in ordered:
        cumulative+=length
        if cumulative>=halfway: n50=length; break
    result={"schema_version":"6.0","reads":len(lengths),"total_bases":sum(lengths),"n50":n50,"minimum_length":min(lengths) if lengths else None,"maximum_length":max(lengths) if lengths else None,"status":"COMPLETE" if lengths else "COMPLETE_EMPTY_RESULT"}
    args.output.write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")


def cmd_depth_manifest(args) -> None:
    sources={row["dataset_id"]:row for row in read_tsv(args.source_manifest)}
    records=[]
    for dataset in args.dataset_id:
        if dataset not in sources: raise ValueError(f"unknown source dataset: {dataset}")
        for fraction in sorted(set(args.fractions)):
            if not 0 < fraction <= 1: raise ValueError("depth fractions must be in (0,1]")
            label=f"{fraction:.2f}".rstrip("0").rstrip(".").replace(".","p")
            records.append({"derived_id":f"{dataset}_d{label}","source_dataset_id":dataset,"population":"full source FASTQ","fraction":f"{fraction:.2f}","seed":args.seed,"method":"SHA-256 read-identifier threshold","biological_replicate_claim":"false","status":args.status})
    write_tsv(args.output,["derived_id","source_dataset_id","population","fraction","seed","method","biological_replicate_claim","status"],records)


def cmd_threshold(args) -> None:
    rows = read_tsv(args.input); outputs = []
    for threshold in args.thresholds:
        retained = []
        unavailable = 0
        for row in rows:
            try: value = float(row[args.support_column])
            except (KeyError, ValueError): unavailable += 1; continue
            if value >= threshold: retained.append(row)
        destination = args.output_dir / f"support_ge_{threshold:g}.tsv"
        write_tsv(destination, list(rows[0].keys()) if rows else [], retained)
        outputs.append({"threshold": threshold, "input_models": len(rows), "retained_models": len(retained),
                        "missing_support": unavailable, "support_semantics": args.support_semantics,
                        "status": "COMPLETE" if unavailable == 0 else "COMPLETE_WITH_UNAVAILABLE_SUPPORT"})
    write_tsv(args.summary, ["threshold", "input_models", "retained_models", "missing_support", "support_semantics", "status"], outputs)


def cmd_support_benchmark(args) -> None:
    predicted=models_from_standardized(args.input); truth=models_from_gtf(args.truth)
    if args.truth_contig_prefix:
        predicted=[model for model in predicted if model["chrom"].startswith(args.truth_contig_prefix)]
        truth=[model for model in truth if model["chrom"].startswith(args.truth_contig_prefix)]
    records=[]
    for threshold in args.thresholds:
        retained=[]; missing=0
        for model in predicted:
            try: support=float(model[args.support_column])
            except (KeyError,ValueError): missing+=1; continue
            if support>=threshold: retained.append(model)
        metrics,_matches=score_truth(retained,truth,args.end_tolerance)
        for level in ("exact_all","exact_multi_exon","exact_single_exon","intron_chain","splice_junction"):
            records.append({"threshold":threshold,"level":level,**metrics[level],"retained_models":len(retained),"missing_support":missing,"support_semantics":args.support_semantics})
    write_tsv(args.output,["threshold","level","tp","fp","fn","precision","recall","f1","status","retained_models","missing_support","support_semantics"],records)


def jaccard(left: set, right: set) -> dict:
    union = left | right
    return {"intersection": len(left & right), "union": len(union), "jaccard": len(left & right) / len(union) if union else None,
            "status": "BOTH_EMPTY" if not union else "COMPLETE"}


def cmd_replicates(args) -> None:
    sets = {}
    for item in args.structure:
        replicate, path = item.split("=", 1)
        models = models_from_standardized(Path(path))
        sets[replicate] = {
            "exact": {(m["chrom"], m["strand"], m["introns"], m["start"], m["end"]) for m in models},
            "intron_chain": {(m["chrom"], m["strand"], m["introns"]) for m in models if m["introns"]},
            "splice_junction": {(m["chrom"], m["strand"], *junction) for m in models for junction in m["introns"]},
        }
    rows = []
    names = sorted(sets)
    for index, left in enumerate(names):
        for right in names[index + 1:]:
            for level in ("exact", "intron_chain", "splice_junction"):
                rows.append({"replicate_left": left, "replicate_right": right, "level": level, **jaccard(sets[left][level], sets[right][level])})
    membership = defaultdict(int)
    for structure in set().union(*(sets[name]["exact"] for name in names)) if names else set():
        membership[sum(structure in sets[name]["exact"] for name in names)] += 1
    write_tsv(args.output, ["replicate_left", "replicate_right", "level", "intersection", "union", "jaccard", "status"], rows)
    args.membership.write_text(json.dumps({"replicates": names, "exact_structures_by_replicate_support": {str(k): v for k, v in sorted(membership.items())}}, indent=2) + "\n", encoding="utf-8")


def ranks(values: list[float]) -> list[float]:
    order = sorted(range(len(values)), key=values.__getitem__); result = [0.0] * len(values); index = 0
    while index < len(order):
        end = index + 1
        while end < len(order) and values[order[end]] == values[order[index]]: end += 1
        rank = (index + end - 1) / 2 + 1
        for position in order[index:end]: result[position] = rank
        index = end
    return result


def correlation(left: list[float], right: list[float]):
    if len(left) < 2: return None
    left_mean, right_mean = statistics.fmean(left), statistics.fmean(right)
    numerator = sum((a-left_mean)*(b-right_mean) for a,b in zip(left,right))
    denominator = math.sqrt(sum((a-left_mean)**2 for a in left)*sum((b-right_mean)**2 for b in right))
    return numerator / denominator if denominator else None


def cmd_quantification(args) -> None:
    truth = {row[args.id_column]: float(row[args.truth_column]) for row in read_tsv(args.truth)}
    predicted = {row[args.id_column]: float(row[args.predicted_column]) for row in read_tsv(args.predicted)}
    ids = sorted(set(truth) | set(predicted)); t = [truth.get(i, 0.0) for i in ids]; p = [predicted.get(i, 0.0) for i in ids]
    rows = [{"transcript_id": identifier, "truth_abundance": a, "predicted_abundance": b,
             "detected": b > 0, "log1p_absolute_error": abs(math.log1p(b)-math.log1p(a)),
             "relative_error": abs(b-a)/a if a else None} for identifier,a,b in zip(ids,t,p)]
    summary = {"schema_version": "6.0", "n_union": len(ids), "n_truth_positive": sum(x > 0 for x in t),
               "n_detected": sum(x > 0 for x in p), "pearson_all_including_zeros": correlation(t,p),
               "spearman_all_including_zeros": correlation(ranks(t), ranks(p)) if ids else None,
               "mean_log1p_absolute_error": statistics.fmean(row["log1p_absolute_error"] for row in rows) if rows else None,
               "mean_relative_error_truth_positive": statistics.fmean(row["relative_error"] for row in rows if row["relative_error"] is not None) if any(row["relative_error"] is not None for row in rows) else None,
               "zero_policy": "union of truth and prediction; missing values represented as zero"}
    write_tsv(args.table, ["transcript_id", "truth_abundance", "predicted_abundance", "detected", "log1p_absolute_error", "relative_error"], rows)
    args.output.write_text(json.dumps(summary, indent=2)+"\n", encoding="utf-8")


def parser() -> argparse.ArgumentParser:
    root = argparse.ArgumentParser(description=__doc__); sub = root.add_subparsers(dest="command", required=True)
    truth = sub.add_parser("truth-match"); truth.add_argument("--prediction", type=Path, required=True); truth.add_argument("--prediction-format", choices=["gtf","standardized"], default="gtf"); truth.add_argument("--truth", type=Path, required=True); truth.add_argument("--truth-contig-prefix",default=""); truth.add_argument("--end-tolerance", type=int, default=0); truth.add_argument("--dataset-id", required=True); truth.add_argument("--caller", required=True); truth.add_argument("--platform", required=True); truth.add_argument("--truth-type", required=True); truth.add_argument("--output", type=Path, required=True); truth.add_argument("--matches", type=Path, required=True); truth.set_defaults(func=cmd_truth)
    inventory=sub.add_parser("truth-inventory"); inventory.add_argument("--truth",type=Path,required=True); inventory.add_argument("--truth-id",required=True); inventory.add_argument("--contig-prefix",default=""); inventory.add_argument("--output",type=Path,required=True); inventory.set_defaults(func=cmd_truth_inventory)
    down = sub.add_parser("downsample"); down.add_argument("--input", type=Path, required=True); down.add_argument("--output", type=Path, required=True); down.add_argument("--fraction", type=float, required=True); down.add_argument("--seed", type=int, required=True); down.add_argument("--manifest", type=Path, required=True); down.set_defaults(func=cmd_downsample)
    stats=sub.add_parser("fastq-stats"); stats.add_argument("--input",type=Path,required=True); stats.add_argument("--output",type=Path,required=True); stats.set_defaults(func=cmd_fastq_stats)
    depth=sub.add_parser("depth-manifest"); depth.add_argument("--source-manifest",type=Path,required=True); depth.add_argument("--dataset-id",action="append",required=True); depth.add_argument("--fractions",nargs="+",type=float,default=[.1,.25,.5,.75,1]); depth.add_argument("--seed",type=int,default=61703); depth.add_argument("--status",default="PLANNED_EXTERNAL_STORAGE"); depth.add_argument("--output",type=Path,required=True); depth.set_defaults(func=cmd_depth_manifest)
    threshold = sub.add_parser("support-threshold"); threshold.add_argument("--input", type=Path, required=True); threshold.add_argument("--support-column", default="support_value"); threshold.add_argument("--support-semantics", required=True); threshold.add_argument("--thresholds", nargs="+", type=float, default=[1,2,3,5,10]); threshold.add_argument("--output-dir", type=Path, required=True); threshold.add_argument("--summary", type=Path, required=True); threshold.set_defaults(func=cmd_threshold)
    support=sub.add_parser("support-benchmark"); support.add_argument("--input",type=Path,required=True); support.add_argument("--truth",type=Path,required=True); support.add_argument("--truth-contig-prefix",default=""); support.add_argument("--support-column",default="support_value"); support.add_argument("--support-semantics",required=True); support.add_argument("--thresholds",nargs="+",type=float,default=[1,2,3,5,10]); support.add_argument("--end-tolerance",type=int,default=0); support.add_argument("--output",type=Path,required=True); support.set_defaults(func=cmd_support_benchmark)
    replicate = sub.add_parser("replicates"); replicate.add_argument("--structure", action="append", required=True, help="REPLICATE=standardized.tsv"); replicate.add_argument("--output", type=Path, required=True); replicate.add_argument("--membership", type=Path, required=True); replicate.set_defaults(func=cmd_replicates)
    quant = sub.add_parser("quantification"); quant.add_argument("--truth", type=Path, required=True); quant.add_argument("--predicted", type=Path, required=True); quant.add_argument("--id-column", default="transcript_id"); quant.add_argument("--truth-column", default="abundance"); quant.add_argument("--predicted-column", default="abundance"); quant.add_argument("--output", type=Path, required=True); quant.add_argument("--table", type=Path, required=True); quant.set_defaults(func=cmd_quantification)
    return root


if __name__ == "__main__":
    arguments = parser().parse_args(); arguments.func(arguments)
