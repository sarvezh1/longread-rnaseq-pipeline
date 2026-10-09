# Input contract

## Sample sheet

The input is a comma-separated file with one row per biological sample and the
following columns:

| Column | Required value | Meaning |
|---|---|---|
| `sample` | Always | Unique identifier using letters, numbers, `.`, `_`, or `-`. |
| `platform` | Always | Exactly `ont` or `pacbio`. |
| `library_type` | Always | One of the implemented types listed below. |
| `preprocessing` | Always | Explicit route; no preprocessing stage is inferred. |
| `condition` | Always | Biological group; tabs and line breaks are prohibited. |
| `replicate` | Always | Replicate identifier; tabs and line breaks are prohibited. |
| `reads` | Always | Existing FASTQ or BAM path appropriate to the route. |
| `primer_fasta` | PacBio HiFi | Explicit primer FASTA for lima/refine; otherwise empty. |
| `primer_pair` | Multiplexed PacBio HiFi | Exact lima output pair, for example `IsoSeqX_bc06_5p--IsoSeqX_3p`; otherwise empty. |
| `ont_kit` | ONT with Pychopper | `PCS109`, `PCS110`, `PCS111`, `PCS114`, `PCB111`, or `PCB114`; otherwise empty. |
| `require_polya` | Column always present | `true` or `false`; meaningful only for PacBio HiFi refine. |

All eleven columns must be present. Fields that do not apply to a route remain
empty. Relative paths are resolved from the directory where Nextflow is
launched. Sample identifiers must be unique.

## Implemented route matrix

| Platform | `library_type` | `preprocessing` | Processing path |
|---|---|---|---|
| `ont` | `cdna` | `pychopper` | NanoPlot → Pychopper → minimap2 splice |
| `ont` | `cdna` | `none` | NanoPlot → minimap2 splice |
| `ont` | `pcr_cdna` | `pychopper` | NanoPlot → Pychopper → minimap2 splice |
| `ont` | `pcr_cdna` | `none` | NanoPlot → minimap2 splice |
| `pacbio` | `isoseq_hifi` | `isoseq` | NanoPlot → lima → refine → cluster2 → pbmm2 |
| `pacbio` | `isoseq_flnc` | `none` | NanoPlot → cluster2 → pbmm2 |

ONT direct RNA, ONT BAM entry, PacBio subreads, and pre-clustered PacBio entry
points are not accepted. Labels outside this matrix fail; they never fall
through to another branch.

## ONT preprocessing choices

Pychopper is applied only when `preprocessing=pychopper`. An explicit supported
`ont_kit` is then required. It orients and trims full-length reads and retains
classified, rescued, unclassified, quality-fail, and length-fail outputs plus
its reports.

Pychopper's tool-native minimum mean quality and fragment-length filters are
explicitly set to zero in the default workflow. This prevents implicit hard
filtering while retaining Pychopper classification. The v1 workflow does not
expose a separate pipeline-level read-length or read-quality filter.

`preprocessing=none` bypasses Pychopper. This is deliberate and recorded in the
normalized sample manifest.

## PacBio preprocessing choices

`isoseq_hifi` expects Q20 HiFi/CCS BAM and does not apply an extra generic
quality filter. `primer_fasta` is required and supplied to both `lima
--isoseq` and `isoseq refine`; primer sequences are never embedded in workflow
code.

One sample-sheet row must resolve through lima to exactly one output BAM. When
a primer FASTA produces multiple barcode-pair BAMs, `primer_pair` must select
the exact pair for that row. If only one BAM is produced, the field may remain
empty. Ambiguous or unmatched output fails explicitly. `--lima_peek_guess true`
exposes lima's optional `--peek-guess` behaviour when appropriate to the
supplied primer set.

`require_polya` controls `isoseq refine --require-polya` and defaults to false
when the cell is empty. The setting is recorded per sample. Cluster2 excludes
singletons by default; `--isoseq_include_singletons true` enables and records
their inclusion for the whole run.

`isoseq_flnc` accepts a BAM already representing FLNC reads. It bypasses lima
and refine but still runs the implemented `cluster2` and pbmm2 stages. A primer
FASTA and `require_polya=true` are invalid for this entry point.

## Reference contract

The v1 target is human GRCh38 with a matching GENCODE annotation. Every run must
provide existing values for:

- `--fasta`;
- `--gtf`;
- `--genome_build`;
- `--annotation_release`.

Nothing is downloaded or inferred. Before preprocessing, the workflow builds a
fresh task-local FASTA index, requires every annotation contig to exist in the
FASTA, checks annotation coordinates against sequence lengths, and records
FASTA, FAI, and annotation SHA-256 checksums. A mismatched pair fails rather
than being silently repaired. Minimap2 and pbmm2 indexes are built once per
active platform branch and reused across samples within a run. The GTF is
required and recorded but does not guide default alignment.

## Annotation-guided alignment

`--annotation_guided_alignment` exists with default `false`. Block 1B does not
implement guided behaviour and rejects `true` with a clear error. This avoids a
partial implementation that guides ONT and PacBio differently and prevents the
annotation from silently biasing novel splice-junction discovery.

## Example

See `assets/samplesheet.example.csv` for all implemented entry points. Its paths
are illustrative and not test data.

## Block 2 input policy

Default mode passes each successfully normalized ONT or PacBio BAM to IsoQuant;
no additional sample-sheet columns are required. Before discovery, Block 2
requires the normalized manifest to match the current sample, platform,
library type, genome build, annotation release, FASTA checksum, and GTF
checksum recorded by Block 1.

`--isoquant_report_novel_unspliced` defaults to null, preserving IsoQuant's
pinned platform-specific policy (false for novel unspliced ONT models and true
for other data types). Set it explicitly to `true` or `false` only when that
scientific choice is intended. The effective policy is recorded per sample.

Block 2 does not add a universal novel-isoform read threshold on top of
IsoQuant's data-type preset.
`--transcriptome_filtering` defaults to false; true is rejected until a
separate, explicit curation policy exists.

Benchmark mode branches from the normalized Block 1 contract to the frozen
Block 3 caller-comparison workflow. Benchmark challengers do not alter the
production transcriptome used by biology mode.

## Block 4 design contract

Every biology-mode sample requires non-empty `sample`, `condition`,
`replicate`, `platform`, and `library_type` metadata. Replicate labels identify
biological replicates within a condition; pseudo-replicates must not be entered.
Platform is metadata and is never inferred to be a biological condition.

Differential analysis is opt-in:

```text
--mode biology
--contrast_numerator treated
--contrast_denominator control
--minimum_biological_replicates 2
--block4_design_formula '~0 + condition'
```

Both contrast labels must exist and differ. The formula must use a no-intercept
parameterization and include `condition`; compatible metadata covariates such
as `platform` or `library_type` may be added only when the resulting design is
full rank. Absent contrasts, missing conditions, insufficient replication,
confounding, and invalid contrasts produce explicit non-run statuses rather
than statistics. The genuine compact ONT and PacBio fixtures have different
biological origins and must never be contrasted with each other.

## Modes and reporting

`nextflow run . --help` prints concise required parameters, modes, common
options, and design safeguards without starting processes. The complete typed
parameter contract is `nextflow_schema.json`.

- `default` runs normalized alignment, production reconstruction, structural
  QC, and final reporting.
- `benchmark` runs the frozen three-caller comparison and final reporting.
- `biology` runs production reconstruction, guarded biological analyses, and
  final reporting.

There is no combined mode in v1 because it would duplicate the frozen
benchmark and production IsoQuant routes. `--report_run_id` may provide a
portable identifier containing letters, numbers, dots, underscores, or
hyphens; otherwise the Nextflow run name is used.

Canonical ONT, PacBio, mixed, and replicated-condition examples are under
`assets/samplesheets/`. Their relative paths are illustrative rather than
bundled data.

## Optional external resources

Large protein homology and domain databases are never downloaded implicitly.
Those optional TransDecoder extensions remain disabled unless an explicit
future configuration supplies the required resources. Human GO enrichment
uses the pinned local annotation package and records its tested background.

## Unsupported in v1

ONT direct RNA, non-human references, single-cell/spatial assays, short-read
reconstruction, de novo transcriptomics, RNA-modification calling, implicit
platform contrasts, pseudo-replicates, automatic sample exclusion, and caller
majority voting are not supported.
