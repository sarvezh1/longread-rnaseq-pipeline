#!/usr/bin/env python3
"""Build deterministic cross-sample matrices from production IsoQuant bundles."""

import argparse
import csv
import hashlib
import json
import re
from collections import defaultdict
from pathlib import Path


def table(path: Path, key: str = "feature_id") -> dict[str, dict]:
    with path.open(newline="") as handle:
        return {row[key]: row for row in csv.DictReader(handle, delimiter="\t")}


def values(path: Path, field: str) -> dict[str, float]:
    return {key: float(row[field]) for key, row in table(path).items()}


def fasta(path: Path) -> dict[str, str]:
    records, name, sequence = {}, None, []
    for line in path.read_text().splitlines():
        if line.startswith(">"):
            if name is not None:
                records[name] = "".join(sequence)
            name, sequence = line[1:].split()[0], []
        elif name is not None:
            sequence.append(line.strip())
    if name is not None:
        records[name] = "".join(sequence)
    return records


def gtf_models(path: Path) -> tuple[dict[str, tuple], dict[str, list[str]]]:
    exons, lines = defaultdict(list), defaultdict(list)
    for line in path.read_text().splitlines():
        if not line or line.startswith("#"):
            continue
        fields = line.split("\t")
        attrs = {}
        for item in fields[8].strip().strip(";").split(";"):
            if item.strip() and " " in item.strip():
                key, value = item.strip().split(" ", 1)
                attrs[key] = value.strip().strip('"')
        tid = attrs.get("transcript_id")
        if tid:
            lines[tid].append(line)
            if fields[2] == "exon":
                exons[tid].append((int(fields[3]), int(fields[4]), fields[0], fields[6]))
    signatures = {tid: tuple(sorted(items)) for tid, items in exons.items()}
    return signatures, lines


def write_matrix(path: Path, features: list[str], samples: list[str], data: dict, integer=False) -> None:
    with path.open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["feature_id", *samples])
        for feature in features:
            row = [data[sample].get(feature, 0.0) for sample in samples]
            writer.writerow([feature, *([int(round(x)) for x in row] if integer else [f"{x:.10g}" for x in row])])


def structural_id(row: dict) -> tuple[str, str]:
    key = "|".join([row.get("chrom", ""), row.get("strand", ""), row.get("exon_coordinates", "")])
    if not all(key.split("|")):
        raise ValueError(f"Transcript {row.get('transcript_id', '<unknown>')} lacks a complete structural identity")
    return "LRT_" + hashlib.sha256(key.encode()).hexdigest()[:16], key


def remap_values(native: dict[str, float], mapping: dict[str, str], sample: str) -> dict[str, float]:
    result = defaultdict(float)
    for transcript_id, value in native.items():
        if transcript_id not in mapping:
            raise ValueError(f"IsoQuant transcript {transcript_id} in {sample} has no Block 2 structural definition")
        result[mapping[transcript_id]] += value
    return dict(result)


def rewrite_gtf_line(line: str, native_id: str, standard_id: str) -> str:
    return re.sub(r'(transcript_id\s+")' + re.escape(native_id) + r'("\s*;)',
                  lambda match: match.group(1) + standard_id + match.group(2), line)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--bundles", required=True, nargs="+", type=Path)
    parser.add_argument("--contrast-numerator", default="")
    parser.add_argument("--contrast-denominator", default="")
    parser.add_argument("--minimum-replicates", type=int, default=2)
    parser.add_argument("--design-formula", default="~0 + condition")
    parser.add_argument("--output", required=True, type=Path)
    args = parser.parse_args()
    args.output.mkdir(parents=True, exist_ok=True)

    bundles = []
    for bundle in args.bundles:
        meta = json.loads((bundle / "sample_metadata.json").read_text())
        bundles.append((meta["sample"], meta, bundle))
    bundles.sort()
    samples = [sample for sample, _, _ in bundles]
    if len(samples) != len(set(samples)):
        raise ValueError("Duplicate samples in Block 4 input")

    gene_counts, transcript_counts, gene_tpm, transcript_tpm = {}, {}, {}, {}
    structures, signatures, selected_gtf, sequences = {}, {}, {}, {}
    source_ids, mapping_rows = defaultdict(set), []
    presence = {}
    for sample, meta, bundle in bundles:
        gene_counts[sample] = values(bundle / "gene_counts.tsv", "count")
        gene_tpm[sample] = values(bundle / "gene_tpm.tsv", "TPM")
        structure = table(bundle / "transcript_structure.tsv", "transcript_id")
        model_signatures, model_lines = gtf_models(bundle / "transcripts.gtf")
        seqs = fasta(bundle / "transcripts.fasta")
        native_to_standard = {}
        for native_id, row in structure.items():
            standard_id, structure_key = structural_id(row); native_to_standard[native_id] = standard_id
            signature = model_signatures.get(native_id, ())
            if standard_id in signatures and signatures[standard_id] != signature:
                raise ValueError(f"Structural ID collision for {standard_id}")
            if standard_id in sequences and native_id in seqs and sequences[standard_id] != seqs[native_id]:
                raise ValueError(f"Structurally identical transcript {standard_id} has incompatible sequences")
            if standard_id in structures:
                prior_gene=structures[standard_id].get("gene_id", "")
                new_gene=row.get("gene_id", "")
                if prior_gene not in {"", "NA"} and new_gene not in {"", "NA"} and prior_gene != new_gene:
                    raise ValueError(f"Structurally identical transcript {standard_id} has incompatible gene assignments")
            standardized=dict(row); standardized["transcript_id"]=standard_id
            structures.setdefault(standard_id, standardized)
            signatures.setdefault(standard_id, signature)
            selected_gtf.setdefault(standard_id, [rewrite_gtf_line(line, native_id, standard_id) for line in model_lines.get(native_id, [])])
            if native_id in seqs: sequences.setdefault(standard_id, seqs[native_id])
            source_ids[standard_id].add(native_id)
            mapping_rows.append({"sample":sample,"caller_transcript_id":native_id,"standard_transcript_id":standard_id,"structure_key":structure_key})
        transcript_counts[sample] = remap_values(values(bundle / "transcript_counts.tsv", "count"), native_to_standard, sample)
        transcript_tpm[sample] = remap_values(values(bundle / "transcript_tpm.tsv", "TPM"), native_to_standard, sample)
        presence[sample] = set(transcript_counts[sample])

    genes = sorted(set().union(*(set(x) for x in gene_counts.values())))
    transcripts = sorted(set().union(*(set(x) for x in transcript_counts.values())))
    write_matrix(args.output / "gene_counts.tsv", genes, samples, gene_counts, True)
    write_matrix(args.output / "transcript_counts.tsv", transcripts, samples, transcript_counts, True)
    write_matrix(args.output / "gene_tpm.tsv", genes, samples, gene_tpm)
    write_matrix(args.output / "transcript_tpm.tsv", transcripts, samples, transcript_tpm)

    with (args.output / "sample_metadata.tsv").open("w", newline="") as handle:
        fields = ["sample", "condition", "replicate", "platform", "library_type"]
        writer = csv.DictWriter(handle, fields, delimiter="\t", lineterminator="\n")
        writer.writeheader(); writer.writerows([{key: meta[key] for key in fields} for _, meta, _ in bundles])
    with (args.output / "transcript_presence.tsv").open("w", newline="") as handle:
        writer = csv.writer(handle, delimiter="\t", lineterminator="\n")
        writer.writerow(["transcript_id", *samples])
        for tid in transcripts:
            writer.writerow([tid, *[int(tid in presence[sample]) for sample in samples]])
    with (args.output / "transcript_id_map.tsv").open("w", newline="") as handle:
        fields=["sample","caller_transcript_id","standard_transcript_id","structure_key"]
        writer=csv.DictWriter(handle,fields,delimiter="\t",lineterminator="\n"); writer.writeheader(); writer.writerows(sorted(mapping_rows,key=lambda row:(row["sample"],row["caller_transcript_id"])))

    metadata_fields = ["transcript_id", "gene_id", "chrom", "strand", "transcript_start", "transcript_end",
                       "transcript_length", "exon_count", "exon_class", "exon_coordinates", "intron_chain",
                       "structural_category", "associated_gene", "associated_transcript", "known_novel", "caller_transcript_ids"]
    with (args.output / "transcript_metadata.tsv").open("w", newline="") as handle:
        writer = csv.DictWriter(handle, metadata_fields, delimiter="\t", lineterminator="\n", extrasaction="ignore")
        writer.writeheader()
        for tid in transcripts:
            row = dict(structures.get(tid, {})); row["transcript_id"] = tid
            row["known_novel"] = "known" if row.get("structural_category") == "full-splice_match" else "novel_or_other"
            row["caller_transcript_ids"] = ",".join(sorted(source_ids[tid]))
            writer.writerow(row)

    with (args.output / "production_transcripts.gtf").open("w") as handle:
        handle.write("# caller-independent Block 4 input: production IsoQuant models only\n")
        for tid in transcripts:
            for line in selected_gtf.get(tid, []): handle.write(line + "\n")
    with (args.output / "production_transcripts.fasta").open("w") as handle:
        for tid in transcripts:
            if tid in sequences: handle.write(f">{tid}\n{sequences[tid]}\n")

    numerator, denominator = args.contrast_numerator, args.contrast_denominator
    conditions = defaultdict(list)
    for _, meta, _ in bundles: conditions[meta["condition"]].append(meta)
    status, reasons = "READY", []
    if not numerator and not denominator:
        status, reasons = "NOT_RUN_NO_EXPLICIT_CONTRAST", ["Both contrast conditions must be supplied explicitly"]
    elif not numerator or not denominator or numerator == denominator:
        status, reasons = "NOT_RUN_INVALID_CONTRAST", ["Contrast numerator and denominator must be distinct and non-empty"]
    elif numerator not in conditions or denominator not in conditions:
        status, reasons = "NOT_RUN_INVALID_CONTRAST", ["A contrast condition is absent from the sample metadata"]
    else:
        for condition in (numerator, denominator):
            replicates = {item["replicate"] for item in conditions[condition]}
            if len(replicates) < args.minimum_replicates:
                reasons.append(f"{condition} has {len(replicates)} biological replicate(s); {args.minimum_replicates} required")
        compared = conditions[numerator] + conditions[denominator]
        for covariate in ("platform", "library_type"):
            mapping = defaultdict(set)
            for item in compared: mapping[item["condition"]].add(item[covariate])
            if len({item[covariate] for item in compared}) > 1 and all(len(x) == 1 for x in mapping.values()) and len({next(iter(x)) for x in mapping.values()}) == 2:
                reasons.append(f"condition is perfectly confounded with {covariate}")
        if reasons: status = "NOT_RUN_INSUFFICIENT_REPLICATION" if any("replicate" in x for x in reasons) else "NOT_RUN_CONFOUNDED_DESIGN"
    design = {"status": status, "contrast": {"numerator": numerator or None, "denominator": denominator or None},
              "minimum_biological_replicates": args.minimum_replicates, "reasons": reasons,
              "design_formula": args.design_formula, "platform_is_biological_condition": False}
    (args.output / "design_status.json").write_text(json.dumps(design, indent=2) + "\n")


if __name__ == "__main__":
    main()
