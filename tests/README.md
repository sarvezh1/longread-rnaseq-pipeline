# Tests

Deterministic structural fixtures and stub contracts live under
`tests/data/structural/` and `tests/contracts/`. They test metadata validation,
platform routing, bypass behaviour, normalized and transcriptome contracts,
reference identity, GTF structure validation, single-exon preservation,
SQANTI3 category summaries, and workflow composition without invoking
scientific containers.

`tests/data/benchmark/` and `tests/contracts/run_block3_contracts.sh` add
clearly synthetic structures for exact, intron-chain-only, partial-junction,
single-exon, caller-specific, and empty-set cases. They are logic fixtures, not
biological validation data.

Real scientific-tool validation lives under `tests/integration/`. Its
preparation script fetches bounded portions of documented public ONT and
PacBio sources, creates deterministic compact subsets, and records SHA-256
checksums. Downloaded reads and runtime outputs remain local and ignored.

The real layer validates Block 1/2 and benchmark mode on the same normalized
BAMs. `run_benchmark.sh` and `check_benchmark_outputs.py` cover all three
callers, common SQANTI3, parseable standardized structures, sparse/empty
contracts, provenance, and published paths. It does not assess biological
accuracy or caller performance.
Block 4 adds `contracts/run_block4_contracts.sh`. It combines structural-ID and
status unit checks, deliberately labelled deterministic statistical fixtures,
pinned-container execution for edgeR, satuRn, SUPPA, TransDecoder, and
IsoformSwitchAnalyzeR, plus a complete biology-mode stub route. The synthetic
counts and transcript structures test plumbing only; they are not biological
validation data.

Block 5 adds deterministic final-report regressions and the read-only
release-readiness checker. The maintained hierarchy is:

1. Level 1: Docker-free schema, parser, contract, report, local-link, and
   repository-hygiene checks using labelled synthetic fixtures.
2. Level 2: Nextflow stub and pinned-container module contracts using minimal
   deterministic data.
3. Level 3: manually triggered genuine compact ONT/PacBio integration. These
   fixtures are provenance documented and are not downloaded by routine CI.

Block 6 adds `contracts/run_block6_validation_tests.py`. It uses deliberately
labelled deterministic structures/counts to test truth matching, TP/FP/FN and
undefined metrics, coordinate tolerance, nested downsampling, support filters,
replicate membership, quantification zeros, provenance, resource parsing, and
publication-table/figure generation. It does not substitute for the deferred
full LRGASP execution documented in `docs/validation.md`. Its `reads.fastq`
fixture consists of eight named synthetic four-base records (32 total bases)
used only to verify deterministic nested downsampling and FASTQ statistics.
