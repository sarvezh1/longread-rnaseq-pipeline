# Block 6 validation framework

This directory is separate from the frozen production workflow. It scores
standardized caller outputs, creates deterministic technical subsets, records
complete benchmark provenance, and renders a standalone offline validation
report. It never changes production defaults.

The source manifest selects matched WTC11 cDNA triplicates: ENCODE experiment
`ENCSR539ZXJ` (ONT) and `ENCSR507JOF` (PacBio), using biosamples
`ENCBS944CBA`, `ENCBS593PKA`, and `ENCBS474NOC`. Because the platforms use
different library and sequencing protocols, comparisons are described as
platform/protocol effects. The supplied PacBio FASTQ files are documented FLNC
reads; a full external run must record the validation-only alignment adapter if
they are not entering through the production PacBio BAM contract.

## Feasibility gate

Run before fetching data:

```bash
python bin/audit_validation_storage.py \
  --manifest validation/manifests/source_datasets.tsv \
  --output validation/resources/storage_audit.json \
  --target /external/benchmark
```

The selected six FASTQs total 55,589,959,267 compressed bytes. The recorded
local audit found less free space than the compressed sources alone and
therefore marked full execution `FULL_LOCAL_EXECUTION_DEFERRED_STORAGE`.
No bounded prefix is represented as a full benchmark.

## Reproducible full-scale route

On external storage, fetch one manifest record at a time, calculate and retain
the generated SHA-256 record, run the frozen production/benchmark processes,
then score their standardized structures:

```bash
python bin/fetch_validation_data.py \
  --manifest validation/manifests/source_datasets.tsv \
  --dataset-id wtc11_ont_cdna_rep1 \
  --output-dir /external/benchmark/source

nextflow run validation/workflow.nf \
  -c validation/nextflow.config \
  --input /external/benchmark/scoring.tsv \
  --outdir /external/benchmark/results \
  --end_tolerance 0 -resume
```

The `scoring.tsv` contract is shown in `manifests/scoring.example.tsv`. Hidden
SIRV truth must be used only by the scoring workflow, never as reconstruction
annotation. For intentionally withheld transcript experiments, the
reconstruction annotation and evaluation truth must be separate files.

Depth subsets use one fixed seed and a SHA-256 threshold on read identifiers,
which makes the 10%, 25%, 50%, 75%, and 100% populations nested. They are
technical subsets, not biological replicates:

```bash
python bin/validation_metrics.py downsample \
  --input source.fastq.gz --output depth_25.fastq.gz \
  --fraction 0.25 --seed 61703 --manifest depth_25.json
```

Support thresholds are evaluated using caller-native semantics. Thresholds are
not equated across callers unless the support definitions are demonstrably
compatible. Use `support-benchmark` to emit TP/FP/FN, precision, recall, and F1
for each predefined threshold.

The offline report is regenerated solely from machine-readable manifests and
metric files:

```bash
python bin/build_validation_report.py \
  --sources validation/manifests/source_datasets.tsv \
  --truth validation/manifests/truth_resources.tsv \
  --runs validation/manifests/benchmark_runs.tsv \
  --storage-audit validation/resources/storage_audit.json \
  --metrics-dir /external/benchmark/results/sirv \
  --output-dir validation
```

Raw and derived reads belong under ignored `validation/data`, `validation/work`,
or external storage. Only manifests, checksums, source tables, scripts, figures,
summary tables, and the report are release artifacts.
