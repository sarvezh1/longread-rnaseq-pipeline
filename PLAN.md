# Development plan

This plan separates production pipeline functionality from optional experiments
that evaluate scientific behaviour. Block 1 is complete and frozen. Block 2
implements the default reconstruction/classification layer. Block 3 implements
the comparative caller benchmark. All three blocks are complete and frozen.

## Block 1 - Inputs, QC, preprocessing, and alignment (complete)

### Block 1A - Bootstrap

- Define the ONT/PacBio sample-sheet and explicit reference contracts.
- Establish modular DSL2 layout, configuration, documentation, and test
  scaffolding.
- Keep default and benchmark concepts separate.

### Block 1B - Platform preprocessing and alignment

- Run ONT raw QC, optional Pychopper, minimap2 splice alignment, sorting,
  indexing, and shared alignment QC.
- Run PacBio raw QC, lima, refine, cluster2, pbmm2 ISOSEQ alignment, indexing,
  and shared alignment QC; support a post-refine FLNC entry point.
- Emit normalized aligned-read records with reference and tool provenance.
- Keep default alignment annotation-unguided.

### Block 1C - Real execution validation

- Prepare local, reproducible genuine ONT PCR-cDNA and PacBio HiFi fixtures.
- Validate all pinned tools, BAM/BAI integrity, QC, reference compatibility,
  checksums, normalized manifests, and the FLNC bypass semantics.

Remaining Block 1 ideas are not active behaviour: a symmetric optional
annotation-guided route and an audited pre/post-QC read-filtering subworkflow.

## Block 2 - Isoform reconstruction and classification (implemented)

- Consume normalized Block 1 BAM/BAI records without realignment.
- Validate Block 1 contract completeness and exact reference identity.
- Run IsoQuant 4.0.0 as the primary per-sample bulk reconstruction engine with
  platform-aware ONT/PacBio CCS presets.
- Retain transcript models, read/gene/transcript assignments, discovered and
  reference counts/TPM, splice-junction/exon information, summaries, and
  runtime provenance.
- Run SQANTI3 6.0.2 as an independent structural classification/QC stage.
- Preserve native SQANTI3 6.x categories and structural QC artifacts.
- Emit raw GTF, transcript FASTA, structural-identity table, machine-readable
  and human-readable summaries, and a standardized transcriptome manifest.
- Preserve single-exon models reported by the tools and document IsoQuant's
  native ONT special treatment.
- Keep raw discovery unfiltered; no curated transcriptome is emitted until an
  explicit curation policy is implemented.
- Preserve per-sample outputs and provide a reversible future hook for unified
  transcriptomes and cross-sample/condition support.

Block 2 validation uses structural/stub tests and the genuine compact Block 1
fixtures. Compact results demonstrate interoperability and output integrity,
not biological sensitivity or accuracy.

## Block 3 - Comparative benchmark mode (implemented)

- Reuse Block 2 IsoQuant and common SQANTI3 modules; add FLAIR 3.0.1 and Bambu
  3.14.0 on the same normalized per-sample BAM/BAI and references.
- Preserve caller-native GTF/BED/FASTA/read-map or abundance/support artifacts
  where generated, plus exact runtime provenance.
- Normalize transcript structures independently of caller IDs and keep
  single-exon logic separate from multi-exon intron-chain logic.
- Report exact structure, intron-chain, individual-junction, transcript-end,
  and gene-assignment concordance as distinct levels.
- Emit three-way membership, SQANTI3-stratified overlap, raw intersection and
  union counts, Jaccard values with explicit empty-set handling, and
  strand-aware configurable TSS/TES bins.
- Preserve sample, condition, replicate, and ONT/PacBio platform metadata for
  future aggregate experiments.
- Do not create a consensus transcriptome or rank callers.

Benchmark mode is an explicit scientific comparison, not a more expensive
synonym for default mode. Optional platform-native challengers remain a future
extension; they are not part of the three-caller Block 3 contract.

## Block 4 - Biological downstream analysis (implemented)

- Build raw-count and IsoQuant-abundance matrices from production models only;
  harmonize run-local transcript labels by exact structure and preserve the map.
- Produce descriptive expression QC without automatic sample exclusion.
- Run edgeR quasi-likelihood gene/transcript expression and satuRn DTU only for
  explicit, replicated, non-confounded designs; preserve filtering and status.
- Generate SUPPA release 2.4 events, PSI, and guarded differential splicing.
- Summarize strand-aware observed starts/ends and configurable putative TSS/TES
  clusters without claiming orthogonal validation.
- Run TransDecoder 6.0.0 with separate `LongOrfs` and `Predict` stages. Retain
  candidate ORFs independently, distinguish `NO_ORF_CANDIDATE` from
  `ORF_CANDIDATE_AVAILABLE`, and record sparse-input prediction unavailability
  without fabricating a primary ORF or asserting non-coding status.
- Run IsoformSwitchAnalyzeR/satuRn switches and keep unsupported functional
  consequences unevaluated; domain/homology hooks remain optional.
- Provide optional local human GO enrichment with the tested gene universe and
  preserved identifier-mapping failures.
- Emit an integrated evidence table and machine-readable Block 4 contract that
  separate observed, statistical, computational, and validated evidence.

Real compact ONT/PacBio fixtures are descriptive execution tests only and are
never contrasted. RNA modification calling remains outside v1.

## Block 5 - Reporting and reproducibility (implemented)

- Generate deterministic offline HTML plus JSON/TSV exports from standardized
  Blocks 1–4 contracts.
- Aggregate supported standard QC with pinned MultiQC 1.35 and keep automated
  summaries disabled.
- Consolidate observed versions, container identities, parameters, references,
  sample definitions, checksums, statuses, and an archive-ready manifest.
- Enforce a common status vocabulary and required-versus-optional completion
  policy without treating valid sparse, empty, or not-run outcomes as failure.
- Provide schema/configuration cleanup, example sheets, CI scaffolding, three
  test levels, release-readiness checks, and public documentation.

## Block 6 - Scientific validation framework (implemented; full data deferred)

- Freeze truth, exact/intron-chain/junction/end definitions before scoring.
- Select matched LRGASP WTC11 ONT/PacBio cDNA triplicates and authoritative
  SIRV-Set 4 truth with explicit accessions, checksums, protocols, and roles.
- Provide deterministic nested depth sampling, caller-native support-threshold
  evaluation, replicate membership/Jaccard, abundance metrics, resource tables,
  publication-source SVGs, and a standalone offline validation report.
- Keep reconstruction annotation separate from hidden evaluation truth and keep
  all Block 1-5 production defaults frozen.
- Record H1-mix quantification, official simulation, and genuine replicated
  H1-hES versus endodermal differential routes without manufacturing results.

The required six WTC11 source FASTQs total 55.59 GB compressed. The audited
local total footprint exceeds available storage, so full caller/platform/depth/
support/replicate/quantification execution is explicitly deferred to the
documented external-storage/HPC route. Framework and deterministic tests are
complete; no unexecuted experiment is reported as an accuracy result.
