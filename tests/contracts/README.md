# Structural and stub contract tests

The path fixtures under `tests/data/structural/` are deliberately empty or
minimal and are suitable only for Nextflow `-stub-run`. In particular,
`ont.fastq` contains one four-base synthetic record solely so the structural
sample sheets resolve to a valid FASTQ. `pacbio_hifi.bam` and
`pacbio_flnc.bam` are valid, minimal BAM files containing one unmapped
four-base synthetic read apiece. They exist only to exercise HiFi preprocessing
and FLNC bypass routing during `-stub-run`; their names and sequences are
invented and they contain no biological data. They were generated from minimal
SAM records with `samtools view --no-PG -b`, which omits variable program
metadata and makes regeneration deterministic. Separate small text fixtures
test the Block 2 GTF, SQANTI3 6.x classification, counts, and contract
utilities. None of these files is biological execution data.

`run_stub_contracts.sh` covers:

1. ONT-only and PacBio-only routes through Block 1 and Block 2;
2. mixed ONT/PacBio channel routing;
3. ONT preprocessing and PacBio FLNC bypasses;
4. IsoQuant-to-SQANTI3-to-standardized-contract workflow composition;
5. unsupported platform/library rejection;
6. missing PacBio primer rejection;
7. missing reference rejection; and
8. transcriptome contract/summary publication.

Before launching Nextflow, it also runs `audit_fixture_paths.py`. The audit
checks every structural sample-sheet path, the fixture groups consumed by the
contract tests, literal CI fixture references, exact `.gitignore` exceptions
for protected read formats, and BAM magic. This makes a missing checkout
fixture fail once with a complete list instead of surfacing route by route.

`run_block2_contracts.sh` covers:

1. Block 1 normalized-contract completeness;
2. genome-build and annotation-release agreement;
3. FASTA/GTF checksum identity;
4. malformed discovered GTF rejection;
5. single-exon transcript preservation;
6. structural coordinates and intron-chain generation; and
7. native SQANTI3 category and support summary parsing.

Run from a POSIX shell with Nextflow and Python 3 available:

```bash
bash tests/contracts/run_block2_contracts.sh
bash tests/contracts/run_stub_contracts.sh
```

Scientific execution validation is kept separately under `tests/integration/`.
`run_block4_contracts.sh` validates abundance-matrix structural harmonization,
explicit non-run states, expression/usage/splicing mathematics, strand-aware
ends, ORF status semantics, and Block 4 DSL2 routing. It uses the locally built
pinned Block 4 containers in addition to Python/Nextflow checks.

`run_block5_report_tests.py` covers complete and zero-transcript samples,
all-empty caller benchmarks, sparse ORF prediction, absent contrasts,
insufficient replication, disabled analyses, missing optional artifacts, and a
failed required reference stage. It requires only Python's standard library.
