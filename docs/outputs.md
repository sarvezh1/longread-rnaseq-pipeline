# Outputs

## Layout

```text
<outdir>/
|-- qc/
|   |-- raw/<sample>.nanoplot/
|   |-- alignment/<sample>.alignment_qc/
|   `-- transcript_structure/<sample>.*
|-- preprocessing/{ont,pacbio}/
|-- alignments/{index,ont,pacbio}/
|-- isoforms/
|   |-- isoquant/per_sample/<sample>.isoquant/
|   `-- sqanti3/per_sample/<sample>.sqanti3/
|-- transcriptome/
|   `-- raw/<sample>.raw.{gtf,fasta}
|-- benchmark/
|   |-- {isoquant,flair,bambu}/
|   |-- sqanti3/{isoquant,flair,bambu}/
|   |-- concordance/{transcripts,intron_chains,splice_junctions,transcript_ends,membership}/
|   `-- summary/
|-- biology/
|   |-- abundance/
|   |-- expression_qc/
|   |-- differential_expression/
|   |-- transcript_usage/
|   |-- alternative_splicing/{events,psi,differential,native}/
|   |-- transcript_ends/
|   |-- orfs/
|   |-- isoform_switches/
|   |-- enrichment/
|   |-- integrated/
|   `-- biology_contract.json
|-- provenance/
|   |-- reference/
|   |-- samples/
|   |-- software_versions/
|   `-- final/
|-- report/
|   |-- final_report.html
|   |-- run_summary.json
|   |-- {samples,analysis_status,tool_versions}.tsv
|   `-- multiqc/
`-- pipeline_info/
```

Block 1 raw-read, preprocessing, alignment, alignment-QC, and normalized
manifest outputs retain their validated meanings. Large FASTQ/BAM outputs use
`--large_file_publish_mode link` by default.

## IsoQuant

`isoforms/isoquant/per_sample/<sample>.isoquant/` retains the complete relevant
IsoQuant result directory, including where generated:

- `*.transcript_models.gtf`: raw discovered known and novel models;
- `*.extended_annotation.gtf`: reference annotation plus novel models;
- `*.transcript_model_reads.tsv.gz`: reads contributing to discovered models;
- `*.read_info.tsv.gz`: read-level gene/transcript assignments;
- discovered transcript/gene counts and TPM tables;
- reference transcript/gene counts and TPM tables;
- exon, splice-junction, and intron-retention count tables;
- SQANTI-like comparison output and execution summary.

The SQANTI-like IsoQuant table is retained but does not replace independent
SQANTI3 classification.

## SQANTI3

`isoforms/sqanti3/per_sample/<sample>.sqanti3/` contains the pinned SQANTI3
6.0.2 classification, junction, corrected GTF/FASTA, parameters, logs, and
structural intermediate files produced by the tool. The optional SQANTI3 R
HTML report is skipped because version 6.0.2 fails on valid zero-junction
classifications; the standardized JSON and Markdown structural summaries are
always generated downstream. The classification table preserves exact SQANTI3
6.x column names and category terms.

When there are no discovered transcript records, the directory instead records
`NO_TRANSCRIPT_MODELS` and contains schema-bearing empty tables. SQANTI3 is not
asked to classify an empty transcriptome.

## Standardized transcriptome

`transcriptome/raw/<sample>.raw.gtf` is the unfiltered primary IsoQuant GTF.
`<sample>.raw.fasta` is the SQANTI3 genome-corrected transcript FASTA when
models exist. No curated output is currently produced.

`qc/transcript_structure/` contains:

- `*.transcript_structure.tsv`: chromosome, strand, exon coordinates, intron
  chain, starts/ends, length, exon count/type, native SQANTI3 category,
  reference assignment, and available support count;
- `*.transcript_summary.json`: machine-readable category and distribution
  summary;
- `*.transcript_summary.md`: concise human-readable summary.

Summary category keys are the native SQANTI3 terms. FSM/ISM and NIC/NNC
aggregates are supplemental counts, not replacements for the native labels.

## Provenance and contracts

`provenance/samples/<sample>.block2_input_validation.json` records successful
Block 1 contract/reference identity checks.

`provenance/samples/<sample>.transcriptome_contract.json` is the authoritative
Block 2 manifest. It records:

- sample, platform, library type, condition, and replicate;
- genome build, annotation release, and reference checksums;
- IsoQuant data type, bulk mode, full-length flag, single-exon policy, and lack
  of pipeline-added filtering;
- raw GTF, transcript FASTA, classification, junction, support, counts,
  structural table, and summary paths;
- links to the Block 1 manifest and observed IsoQuant/SQANTI3 version records;
- the structural definition of transcript identity.

The runtime channel exposes the tuple documented in
[architecture.md](architecture.md). Paths in the JSON are relative to the run
output directory.

## Benchmark outputs

`benchmark/isoquant`, `benchmark/flair`, and `benchmark/bambu` retain native
caller results. FLAIR includes GTF, BED, transcript FASTA, read-to-isoform map,
support table, and execution log when models exist; legitimate no-junction
inputs retain empty versions of those artifacts and an explicit status. Bambu
retains its discovery R object, all-model and read-supported GTFs,
caller-native novelty/read-support table, abundance when available,
quantification status, and execution log.

`benchmark/sqanti3/<caller>/` holds identically configured SQANTI3 results for
each caller. The corrected FASTA is present even for empty outputs (zero-byte
sequence content with schema-bearing classification/junction tables).

Standardized `*.standardized_transcripts.tsv` tables contain sample metadata,
caller and caller ID, chromosome/strand, exons, intron chain, genomic ends,
strand-aware TSS/TES, length, gene assignment, SQANTI3 category, and native
support value/semantics. Junction tables preserve one row per model junction;
canonical/reference/support fields remain `NA` where the caller/common output
does not provide a genuinely comparable value.

Per-sample concordance outputs include:

- `*.membership.tsv`: structural groups and all-three/pairwise/caller-only membership;
- `*.pairwise_metrics.tsv`: intersection, union, Jaccard, and empty-set status for exact multi-exon structures, intron chains, junctions, and exact single-exon coordinates;
- `*.transcript_ends.tsv`: strand-aware signed/absolute TSS/TES differences and bins for shared intron chains;
- `*.partial_junction_overlap.tsv`: transcript pairs sharing only part of a junction chain;
- `*.gene_overlap.tsv`: Level D common gene assignments;
- `*.benchmark_summary.json` and `*.benchmark_contract.json`: caller counts, categories, memberships, metrics, parameters, platform metadata, published paths, and runtime provenance.

## Block 4 biology outputs

`biology/abundance/` is the analysis-ready source of truth: raw gene and
transcript count matrices, IsoQuant gene/transcript TPM matrices, sample
metadata, design status, production GTF/FASTA, SQANTI3-aware transcript
metadata, per-sample presence, and the reversible caller-ID to structural-ID
map. Statistical methods never substitute TPM for counts.

`expression_qc/` contains library totals, detected-feature counts,
sample correlations, SVG summaries, and PCA coordinates when sample/feature
dimensions permit. It does not automatically exclude samples.

`differential_expression/` keeps gene and transcript edgeR tables, tested/not
tested filtering tables, method status, and versions. Transcript results retain
gene, SQANTI3, and known/novel context. `transcript_usage/` separately contains
per-sample within-gene usage, satuRn results, and explicit reasons why features
were not testable.

`alternative_splicing/` contains SUPPA-native event definitions, normalized PSI
with `NA` preserved, differential ΔPSI/p-values/BH-adjusted values when the
design is valid, and native execution artifacts. Expression, usage, and PSI are
not interchangeable.

`transcript_ends/` records genomic ends, strand-aware putative TSS/TES,
unmodified coordinates, configurable cluster IDs, support, abundance, and
per-gene summaries. These are long-read observations, not validated promoter
or polyadenylation sites.

`orfs/` contains all selected and alternative TransDecoder candidates,
coordinates, orientation, completeness, peptide lengths/sequences, and explicit
per-transcript `NO_ORF_CANDIDATE` or `ORF_CANDIDATE_AVAILABLE` states.
`orf_status.json` separately records run-level `PREDICTION_COMPLETE` or
`PREDICTION_NOT_AVAILABLE_SPARSE_INPUT`. In the latter case, `LongOrfs`
candidates are retained, no primary ORF is fabricated, and absence of a
reliable `Predict` result is not labelled as non-coding. `isoform_switches/`
keeps the native analysis object, statistical switch table, and only
evidence-supported consequence fields. Homology/domain claims remain absent
unless optional databases are supplied.

`enrichment/` is schema-bearing and `NOT_RUN_DISABLED` by default. Enabled GO
analysis records tested sets, background universe, and every identifier mapping
failure. `integrated/` joins observed abundance/structure to separate
statistical and computational evidence without producing prose claims.

`biology_contract.json` is the Block 5 interface. It lists all standardized
paths, per-analysis statuses, observed runtime versions, evidence hierarchy,
and the integrated-table checksum.

## Block 5 report and reproducibility outputs

| Output | Format | Required | Evidence role | Empty/not-run behavior |
|---|---|---:|---|---|
| `report/final_report.html` | Offline HTML | Yes | Navigation, QC, and status overview | Shows explicit status; never substitutes zero for not-run |
| `report/run_summary.json` | JSON, schema 5.0 | Yes | Machine-readable run summary | Missing required contracts set overall `FAILED` |
| `report/samples.tsv` | TSV | Yes | Descriptive per-sample QC | A genuine zero-transcript result remains zero |
| `report/analysis_status.tsv` | TSV | Yes | Optional/inferential-stage status | Keeps empty, disabled, unavailable, and invalid design distinct |
| `report/tool_versions.tsv` | TSV | Yes | Expected/observed software audit | Unused mode-specific tools are `NOT_OBSERVED_NOT_RUN` |
| `report/multiqc/multiqc_report.html` | Offline HTML | Yes | Standard sequencing/alignment QC | Remains valid when only sparse supported logs exist |
| `provenance/final/reproducibility_manifest.json` | JSON | Yes | Archive-ready software/reference/parameter definition | Omits ephemeral work paths |
| `provenance/final/provenance_bundle.json` | JSON | Yes | Checksums and completion state | Malformed content fails integrity validation |
| `provenance/final/tool_audit.tsv` | TSV | Yes | Canonical version/container audit | Mismatches and tag-only image pins are explicit |

The custom report consumes standardized contracts. MultiQC is restricted to
formats it supports natively. Neither report embeds BAM, FASTQ, databases, or
external web content.

## Absent outputs

No block produces a majority-vote caller transcriptome, biological truth set,
automatic sample exclusion, experimentally validated promoter/TES/function
claim, or RNA-modification call.

## Block 6 validation outputs

`validation/manifests/` contains source, truth, derived-depth, run, metric, and
deferred differential-design contracts. `validation/truth/` retains only the
small verified truth inventory; downloaded reference/read files remain ignored.
`validation/source_tables/`, `figures/`, and `tables/` are deterministic report
sources. `validation/report/index.html` and `validation_summary.json` distinguish
complete metrics from deferred experiments. Raw FASTQ/BAM, Nextflow work, and
local benchmark results are excluded from the source tree.

Truth tables are evaluative only when tied to an explicit truth resource.
Caller/platform/replicate concordance without truth remains descriptive.
Undefined metrics are null with an explicit status; deferred or not-run analyses
are never rendered as zero effect.
