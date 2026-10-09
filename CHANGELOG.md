# Changelog

The format follows [Keep a Changelog](https://keepachangelog.com/en/1.1.0/)
and the project uses [Semantic Versioning](https://semver.org/).

## [Unreleased]

No changes yet.

## [0.1.0] - 2026-10-08

### Added

- Platform-specific ONT cDNA/PCR-cDNA and PacBio Iso-Seq preprocessing,
  alignment, QC, validation, and normalized BAM/BAI contracts.
- Production IsoQuant 4.0.0 reconstruction and SQANTI3 6.0.2 structural
  characterization with standardized transcriptome contracts.
- Comparative IsoQuant, FLAIR 3.0.1, and Bambu 3.14.0 benchmark mode with
  exact-structure, intron-chain, junction, transcript-end, gene, and membership
  concordance outputs.
- Production biological-analysis mode with abundance matrices, expression QC,
  guarded edgeR DE/DTE, satuRn DTU, SUPPA events/PSI, transcript-boundary
  summaries, TransDecoder ORFs, IsoformSwitchAnalyzeR, and optional GO
  enrichment.
- Explicit statistical and sparse-input status contracts, including separate
  TransDecoder candidate and primary-prediction states.
- MultiQC 1.35 standard-QC aggregation and a deterministic offline final report
  generated from Blocks 1–4 contracts.
- Consolidated tool/container audit, resolved parameters, reference hashes,
  run summary, and archive-ready reproducibility manifest.
- Frozen v1 sample-sheet examples, resource ceilings, local/Docker and
  Apptainer/Singularity profile hooks, CI scaffolding, three-level tests, and a
  read-only release-readiness checker.
- External LRGASP/SIRV validation manifests, a storage feasibility gate,
  deterministic truth/depth/support/replicate/quantification metrics, benchmark
  provenance validation, resource normalization, source-table SVG generation,
  standalone Nextflow scoring, and an offline validation report.

### Scientific policy

- IsoQuant remains the production caller. Benchmark caller agreement is context
  and never creates consensus truth.
- Differential inference requires an explicit estimable replicated design;
  unrelated ONT/PacBio fixtures are never contrasted.
- Transcript starts/ends are putative observations, and ORFs, switching
  consequences, and enrichment remain computational predictions.
- Empty, disabled, invalid-design, and sparse-input outcomes are reported
  explicitly rather than converted to biological zeros or technical failures.

### Validation

- Deterministic structural/statistical fixtures cover routing, concordance,
  expression, usage, splicing, boundaries, ORFs, switching, contracts, final
  reporting, and failure/status edge cases.
- Provenance-documented compact genuine ONT and PacBio fixtures validate local
  Docker execution and output integrity without claiming biological
  performance.
- Authoritative SIRV truth files and checksums were verified locally. Full
  WTC11 and replicated biological runs remain explicitly deferred because
  their audited storage footprint exceeds available local capacity; no
  publication accuracy values are claimed for unexecuted runs.
