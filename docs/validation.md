# Scientific validation and publication benchmarking

## Scope and status

Block 6 is an external validation layer, not a production mode and not a source
of automatic defaults. It distinguishes execution integrity, truth-bearing
reconstruction accuracy, caller concordance, platform/protocol effects, depth,
support thresholds, replicate reproducibility, quantification, real
differential execution, and resource use. No global accuracy score is defined.

The framework and deterministic metric tests are complete. Full external-data
execution is currently `PARTIAL_VALIDATION_DEFERRED_EXTERNAL_STORAGE`: the six
selected WTC11 FASTQs total 55,589,959,267 compressed bytes, while the recorded
workstation audit had 46,352,101,376 bytes free. The estimated unpacked, work,
and result footprint is 611,489,551,937 bytes. This is a feasibility outcome,
not a biological zero and not a completed full-scale benchmark.

## Why LRGASP

The [LRGASP RNA-seq data matrix](https://lrgasp.github.io/lrgasp-submissions/docs/rnaseq-data-matrix.html)
provides documented human long-read datasets, matched biological triplicates,
multiple protocols/platforms, SIRV spike-ins, mixtures, and official simulation
resources. The [LRGASP evaluation documentation](https://lrgasp.github.io/lrgasp-submissions/docs/evaluation.html)
provides the community evaluation context. Source datasets, exact accessions,
URLs, file sizes where verified, roles, and local statuses are frozen in
`validation/manifests/source_datasets.tsv`.

The primary paired design is WTC11 cDNA:

| Platform | Experiment | Replicate FASTQ accessions | Biosamples |
|---|---|---|---|
| ONT | ENCSR539ZXJ | ENCFF263YFG, ENCFF023EXJ, ENCFF961HLO | ENCBS944CBA, ENCBS593PKA, ENCBS474NOC |
| PacBio | ENCSR507JOF | ENCFF563QZR, ENCFF370NFS, ENCFF245IPA | ENCBS944CBA, ENCBS593PKA, ENCBS474NOC |

The pairing is biological, but preparation and sequencing protocols differ.
Results must therefore be called a **platform/protocol effect**, not a pure
platform effect. Replicates are never pooled before reproducibility analysis.

H1-mix cDNA experiments `ENCSR957QYS` (ONT) and `ENCSR731MFY` (PacBio) are
reserved for mixture-aware quantification validation. Official human simulation
records `syn25683375` (ONT-like) and `syn25683376` (PacBio-like) are lower
priority and require the documented Synapse access route.

For real replicated differential execution, the selected deferred design is
ONT cDNA H1-hES (`ENCSR016IKV`) versus endodermal (`ENCSR485VRY`), with three
genuine biological replicates per condition and the same size-selected cDNA
protocol. Its six file accessions and biosamples are in
`validation/manifests/differential_candidate.tsv`. One verified input alone is
12,036,441,352 compressed bytes, making the complete analysis plus work files
infeasible in the current free space. No pseudo-replicates were created.

## Truth resources and leakage protection

The authoritative [LRGASP annotation archive](https://zenodo.org/records/10884185)
supplies SIRV-Set 4 sequence and annotation. Downloaded files are verified
against the published MD5 values. The local inventory confirms 84 SIRV
transcripts: 61 multi-exon, 23 single-exon, 114 unique splice junctions, across
22 SIRV reference sequences. ERCC records present in the combined annotation
are excluded from SIRV transcript scoring by the explicit `SIRV` contig prefix.

Reconstruction reference and evaluation truth are separate concepts. Hidden or
artificial novel truth must never be included in caller annotation and then
scored as novel discovery. GENCODE-relative novelty is not experimental truth.
Where truth is absent, outputs describe novelty and reproducibility, not false
positives.

## Frozen structural definitions

Definitions are machine-readable in `validation/manifests/metric_definitions.json`
and are fixed before comparative scoring.

- Multi-exon exact: chromosome/reference sequence, strand, complete ordered
  intron chain, and strand-aware TSS/TES within one configured tolerance.
- Intron-chain identity: chromosome, strand, and complete ordered chain,
  independent of transcript ends.
- Junction identity: chromosome, strand, intron start, and intron end.
- Single-exon exact: chromosome, strand, and TSS/TES tolerance; single-exon
  models never enter intron-chain metrics.
- End error: signed and absolute strand-aware TSS/TES differences for
  intron-chain-matched models.

The default benchmark tolerance is 0 bp. Matching is deterministic one-to-one
on unique structures. Raw TP, FP, and FN accompany precision, recall, and F1.
Empty truth and zero detection have explicit, distinct statuses and null values
where a metric is mathematically undefined.

## Experiments

Caller experiments use one normalized alignment per dataset/replicate for
IsoQuant, FLAIR, and Bambu. Caller IDs are never identity. Block 3 supplies
exact/intron-chain/junction/end concordance, all-three/pairwise/caller-only
membership, and raw intersection/union with Jaccard. Agreement is not truth.

Depth populations are produced with a fixed SHA-256 threshold over seed 61703
and read identifier at 10%, 25%, 50%, 75%, and 100%. Populations are nested and
are technical subsets, never pseudo-replicates. Each level records reads,
bases, N50, reconstruction/truth metrics, abundance metrics, and runtime.

Support experiments retain caller-native semantics and test thresholds 1, 2,
3, 5, and 10 only where support is available. Equal numeric thresholds are not
declared comparable across callers. Production filtering is unchanged.

Replicate outputs retain pairwise exact-structure, intron-chain, and junction
intersections/unions/Jaccards; transcript detection in 1/3, 2/3, and 3/3;
abundance correlation; SIRV consistency; and transcript-end variability.
Quantification uses the union of truth and prediction and keeps undetected
transcripts as explicit zeros for the documented all-transcript correlations.

Repeated biological summaries should report replicate values, central tendency,
and variability. Paired platform/protocol summaries retain biosample pairing.
Any hypothesis test must declare its hypothesis, effect size, and multiplicity
correction; it cannot establish universal caller superiority.

## Reproduction

`validation/README.md` contains the storage gate, fetch, scoring, and report
commands. Required records are:

- source, truth, derived-data, differential-candidate, and benchmark-run
  manifests under `validation/manifests/`;
- per-run provenance created and validated by `benchmark_provenance.py`;
- deterministic truth/depth/support/replicate/quantification metrics from
  `validation_metrics.py`;
- Nextflow trace normalization from `summarize_validation_resources.py`;
- source-table-driven SVGs from `plot_validation.py`;
- offline report and JSON summary from `build_validation_report.py`.

Every executed benchmark record must include accession/URL/checksum, subset
method and seed, platform/protocol/replicate, reference/annotation hashes,
caller/version/parameters, support/depth, container identity, runtime, and the
metric-definition hash. Records missing these fields cannot enter publication
summaries.

## Limitations

- Full WTC11, H1-mix, real differential, and official simulation executions are
  deferred to external storage/HPC and have no pipeline accuracy results yet.
- The downloaded SIRV files and inventory validate truth provenance and parsing,
  not caller reconstruction performance.
- SIRV spike-ins are controlled constructs and do not capture all properties of
  human transcript complexity.
- Protocol differs with platform in the primary paired design.
- A transcript reproducible across callers or replicates is not necessarily
  biologically correct; a caller-specific or one-replicate model is not
  necessarily false.
- No benchmark result has changed a production default.
