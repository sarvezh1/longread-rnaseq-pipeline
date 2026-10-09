#!/usr/bin/env bash
set -euo pipefail

root_dir=$(cd "$(dirname "${BASH_SOURCE[0]}")/../.." && pwd)
out_dir=${1:-"${root_dir}/.integration-fixtures/results"}
samtools_image='community.wave.seqera.io/library/htslib_samtools:1.24--d697cfb9dce007cd'

test -s "${out_dir}/provenance/reference/reference_validation.tsv"
grep -Fx $'status\tPASS' "${out_dir}/provenance/reference/reference_validation.tsv"
absolute_out=$(cd "${out_dir}" && pwd)

python3 - "${absolute_out}" <<'PY'
import json
import pathlib
import re
import sys

outdir = pathlib.Path(sys.argv[1])
manifests = sorted((outdir / "provenance" / "samples").glob("*.normalized.json"))
if not manifests:
    raise SystemExit("No normalized manifests found")

for manifest_path in manifests:
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for field in ("sorted_bam", "bam_index", "alignment_qc_directory",
                  "runtime_tool_versions", "reference_validation"):
        target = outdir / manifest[field]
        if not target.exists():
            raise SystemExit(f"{manifest_path.name}: missing {field}: {target}")
    versions = (outdir / manifest["runtime_tool_versions"]).read_text(encoding="utf-8")
    if ': ""' in versions or any(line.strip() == "VERSIONS" for line in versions.splitlines()):
        raise SystemExit(f"{manifest_path.name}: malformed or empty runtime version record")

contracts = sorted((outdir / "provenance" / "samples").glob("*.transcriptome_contract.json"))
if len(contracts) != len(manifests):
    raise SystemExit(
        f"Expected one Block 2 contract per normalized manifest; observed {len(contracts)} and {len(manifests)}"
    )

fasta_contigs = set()
reference_fasta = outdir.parent / "reference" / "mini_grch38.fa"
with reference_fasta.open() as handle:
    for line in handle:
        if line.startswith(">"):
            fasta_contigs.add(line[1:].split()[0])

for contract_path in contracts:
    contract = json.loads(contract_path.read_text(encoding="utf-8"))
    if contract.get("contract_version") != "2.0":
        raise SystemExit(f"{contract_path.name}: unexpected contract version")
    outputs = contract.get("outputs", {})
    required = (
        "primary_transcript_annotation", "transcript_fasta", "transcript_classification",
        "splice_junction_classification", "transcript_support", "isoquant_output_directory",
        "sqanti3_output_directory", "structural_identity_table", "structural_summary",
    )
    for field in required:
        target = outdir / outputs.get(field, "")
        if not outputs.get(field) or not target.exists():
            raise SystemExit(f"{contract_path.name}: missing Block 2 {field}: {target}")

    for field in ("block1_normalized_manifest", "isoquant_versions", "sqanti3_versions"):
        target = outdir / contract.get("provenance", {}).get(field, "")
        if not target.exists():
            raise SystemExit(f"{contract_path.name}: missing provenance {field}: {target}")
        if field.endswith("versions"):
            version_lines = target.read_text(encoding="utf-8").splitlines()
            if any(line.strip() and ":" not in line for line in version_lines):
                raise SystemExit(f"{contract_path.name}: malformed version YAML: {target}")

    isoquant_versions = outdir / contract["provenance"]["isoquant_versions"]
    if not re.search(r'^  isoquant: "4\.0\.0"$', isoquant_versions.read_text(), re.MULTILINE):
        raise SystemExit(f"{contract_path.name}: unexpected IsoQuant runtime version")
    sqanti3_versions = outdir / contract["provenance"]["sqanti3_versions"]
    if not re.search(r'^  sqanti3: "6\.0\.2"$', sqanti3_versions.read_text(), re.MULTILINE):
        raise SystemExit(f"{contract_path.name}: unexpected SQANTI3 runtime version")

    classification = outdir / outputs["transcript_classification"]
    header = classification.open().readline().rstrip("\n").split("\t")
    for column in ("isoform", "chrom", "strand", "length", "exons", "structural_category"):
        if column not in header:
            raise SystemExit(f"{contract_path.name}: classification lacks {column}")

    summary = json.loads((outdir / outputs["structural_summary"]).read_text())
    gtf = outdir / outputs["primary_transcript_annotation"]
    transcript_count = 0
    with gtf.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            if not line.strip() or line.startswith("#"):
                continue
            fields = line.rstrip("\n").split("\t")
            if len(fields) != 9:
                raise SystemExit(f"{gtf}: malformed line {line_number}")
            if fields[0] not in fasta_contigs:
                raise SystemExit(f"{gtf}: unknown reference contig {fields[0]}")
            if fields[2] == "transcript":
                transcript_count += 1
    if summary.get("total_transcript_models") != transcript_count:
        raise SystemExit(f"{contract_path.name}: GTF/summary transcript count mismatch")
    transcript_fasta = outdir / outputs["transcript_fasta"]
    if transcript_count > 0 and transcript_fasta.stat().st_size == 0:
        raise SystemExit(f"{contract_path.name}: transcript FASTA is empty despite discovered models")

    fasta_records = {}
    current_id = None
    sequence = []
    with transcript_fasta.open() as handle:
        for line_number, line in enumerate(handle, start=1):
            line = line.strip()
            if not line:
                continue
            if line.startswith(">"):
                if current_id is not None:
                    fasta_records[current_id] = "".join(sequence)
                current_id = line[1:].split()[0]
                if not current_id or current_id in fasta_records:
                    raise SystemExit(f"{transcript_fasta}: invalid or duplicate FASTA ID on line {line_number}")
                sequence = []
            elif current_id is None:
                raise SystemExit(f"{transcript_fasta}: sequence precedes first FASTA header")
            else:
                sequence.append(line)
        if current_id is not None:
            fasta_records[current_id] = "".join(sequence)

    classified_lengths = {}
    with classification.open() as handle:
        columns = handle.readline().rstrip("\n").split("\t")
        isoform_index = columns.index("isoform")
        length_index = columns.index("length")
        for line_number, line in enumerate(handle, start=2):
            if not line.strip():
                continue
            fields = line.rstrip("\n").split("\t")
            classified_lengths[fields[isoform_index]] = int(fields[length_index])

    if set(fasta_records) != set(classified_lengths):
        raise SystemExit(
            f"{contract_path.name}: transcript FASTA/classification ID mismatch: "
            f"{sorted(fasta_records)} != {sorted(classified_lengths)}"
        )
    for isoform, expected_length in classified_lengths.items():
        observed_length = len(fasta_records[isoform])
        if observed_length != expected_length:
            raise SystemExit(
                f"{contract_path.name}: FASTA length for {isoform} is {observed_length}, "
                f"classification reports {expected_length}"
            )
PY

for bam in "${out_dir}"/alignments/ont/*.bam "${out_dir}"/alignments/pacbio/*.bam; do
    [[ -e "${bam}" ]] || continue
    relative_bam=${bam#"${out_dir}"/}
    docker run --rm -v "${absolute_out}:/results:ro" "${samtools_image}" bash -lc "
        set -euo pipefail
        samtools quickcheck -v '/results/${relative_bam}'
        samtools idxstats '/results/${relative_bam}' >/dev/null
        test \"\$(samtools view -c '/results/${relative_bam}')\" -gt 0
        samtools view -H '/results/${relative_bam}' | grep -q '^@HD.*SO:coordinate'
        samtools view -H '/results/${relative_bam}' | grep -q '^@RG.*SM:'
    "
done

printf 'Validated Block 1 BAMs/manifests and Block 2 transcriptome contracts in %s\n' "${out_dir}"
