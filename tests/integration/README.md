# Real-fixture integration tests

This layer executes the pinned scientific containers on genuine long-read RNA-seq
records. It is separate from `tests/contracts/`, which uses stubs and remains the
fast, container-free structural test layer.

## Fixtures

The preparation script creates local, ignored fixtures rather than redistributing
the source reads. Run from a Linux environment with `curl`, `gzip`, `awk`, Docker,
and the pinned samtools image available:

```bash
bash tests/integration/prepare_fixtures.sh
bash tests/integration/run_real.sh
bash tests/integration/check_outputs.sh
docker build -t longread-rnaseq-pipeline/bambu:3.14.0-bioc3.23-r4.6 containers/bambu
bash tests/integration/run_benchmark.sh
```

The ONT fixture is the first 4,000 complete reads in the decompressed 10 MiB
prefix of ENA run SRR13619571 (GEO GSM5033695). The source is human coronary
artery smooth-muscle-cell nuclear RNA prepared with ONT SQK-PCS109 PCR-cDNA.

The PacBio fixture is the first 10,000 complete BAM records recoverable from a
64 MiB byte prefix of PacBio's public UHRR monomer Iso-Seq HiFi BAM. The script
rewrites the records with pinned samtools 1.24 so the result has a complete BAM
EOF block while preserving the source CCS read group and tags. The corresponding
official Iso-Seq v2 12-plex primer FASTA is fetched unchanged.

The compact reference contains the GRCh38.p14 TGFBI interval
`NC_000005.10:136027988-136064818` and the human mitochondrial reference
`NC_012920.1`. GENCODE 49 basic-annotation records for TGFBI and chrM are
extracted, and TGFBI coordinates are shifted onto the local `chr5_TGFBI`
contig. FASTA and GTF sequence names therefore match.

All byte ranges, URLs, transformations, read counts, and checksums are recorded
by `prepare_fixtures.sh` and its generated `SHA256SUMS`. The ENA/SRA and PacBio
source pages do not attach a fixture-specific redistribution licence; therefore
the read fixtures are fetched locally and are not included in this source tree.
Users remain responsible for the source archives' current access and reuse terms.

## Scope

The mixed test exercises ONT PCR-cDNA and PacBio `isoseq_hifi` routing through
raw QC, platform-specific preprocessing, alignment, BAM validation, alignment
QC, and normalized manifests. In default mode it continues through per-sample
IsoQuant bulk discovery, SQANTI3 structural classification when models exist,
structural summaries, and standardized transcriptome manifests.
Benchmark mode reuses the same normalized BAMs and references for IsoQuant,
FLAIR, and Bambu, applies common SQANTI3 characterization, and checks the
published concordance contract.

Biology mode reuses the production IsoQuant/SQANTI3 outputs to validate matrix
construction, descriptive expression QC, transcript-boundary summaries, ORF
execution where models exist, explicit non-run statuses, provenance, and the
standardized Block 4 contract. The ONT and PacBio fixtures have different
biological sources and are never used as conditions in differential inference.
No real replicated differential dataset is bundled; that path is validated by
the deterministic statistical fixture in `tests/data/biology/`.

Block 5 resumes these outputs and checks offline MultiQC, custom HTML/JSON/TSV
reporting, relative links, consolidated provenance, zero-transcript display,
benchmark empty-result display, and the PacBio sparse ORF status. This manual
integration layer is not part of routine CI.

To exercise the `isoseq_flnc` entry point, copy a
genuine `*.flnc.bam` produced by the mixed test's refine task to
`.integration-fixtures/pacbio/pacbio_real.flnc.bam`, then run the generated
`samplesheet.flnc.csv`. FLNC is a semantic contract: no universal BAM tag or
header signature reliably distinguishes every valid third-party post-refine
FLNC BAM from a clustered transcript BAM.

The compact fixtures can yield few or no models under tool-native support
policies. That outcome is recorded and is not worked around by lowering
scientific thresholds. These tests establish execution and file integrity only;
they do not assess biological accuracy, sensitivity, or caller performance.
