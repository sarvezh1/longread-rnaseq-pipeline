# Blocks 1-4 real-execution validation record

Block 1 validation was performed on 2026-10-06 and Block 2 validation on
2026-10-07 with Nextflow 26.04.6 and Docker. These checks establish executable
interoperability and output integrity, not biological accuracy or benchmarking
performance.

## Fixtures and provenance

### Oxford Nanopore

- Source: ENA/SRA run `SRR13619571`, GEO sample `GSM5033695`.
- Biology/library: human coronary artery smooth-muscle-cell nuclear RNA, ONT
  SQK-PCS109 PCR-cDNA, MinION.
- URL: `https://ftp.sra.ebi.ac.uk/vol1/fastq/SRR136/071/SRR13619571/SRR13619571_1.fastq.gz`.
- Full source metadata: 1,320,000 reads, 537,084,062 bytes, ENA MD5
  `e31739d3beae1cbdf3221ae2b8eb4707`.
- Selection: first 4,000 complete reads in the decompressed first 10 MiB byte
  range; sequences and qualities were not modified.
- Fixture SHA-256:
  `0a1cbcc278436d813b6dbfd7d6ce125f8a69c19aff5edf98b04f90700b21fc64`.

### PacBio

- Source: PacBio public Kinnex full-length RNA example, regular non-Kinnex
  12-plex UHRR monomer Iso-Seq HiFi BAM
  `m64307e_230628_025302.hifi_reads.bam`.
- URL: `https://downloads.pacbcloud.com/public/dataset/Kinnex-full-length-RNA/DATA-SQ2-UHRR-Monomer/1-CCS/m64307e_230628_025302.hifi_reads.bam`.
- Source metadata: 4,743,856,041 bytes; server ETag
  `21860edb4a85500c6ed067e4027ccfec-71`.
- Selection: first 10,000 complete BAM records from byte range 0-67,108,863,
  rewritten with samtools 1.24 to add a complete EOF block; CCS read-group and
  record tags were preserved.
- Downloaded-range SHA-256:
  `e8b657564241e0dc37c23ec296bc6c14ab3a67637f1783c437de36c2a993be32`.
- Rewritten-fixture SHA-256:
  `e4bae6eaae5b75b17c1ee56d092fc2496ae1f03898e5e04ce863c31920087a21`.
- Official `IsoSeq_v2_primers_12.fasta` SHA-256:
  `af554b24b45b7cc4184888b3acb7148befa446f746a3f9cc4a9930b94c8e16b4`.

The read-source pages do not state a fixture-specific redistribution licence.
The preparation mechanism fetches reads locally; this project does not
redistribute them.

### Reference

- GRCh38.p14 RefSeq TGFBI interval
  `NC_000005.10:136027988-136064818`, renamed `chr5_TGFBI`.
- Human mitochondrial reference `NC_012920.1`, renamed `chrM`.
- GENCODE 49 basic annotation restricted to TGFBI and chrM; TGFBI coordinates
  were shifted to the local interval.
- Generated FASTA SHA-256:
  `ec093f6c9c5fe5c4293c1c3c8276eccf94002022b1237ae4eaddc2f658f69fc9`.
- Generated GTF SHA-256:
  `af0b5a52498a07347aeffbf05c2b291f62b014390049d9e702736e954d9f4e63`.

The workflow confirmed both GTF contigs were present in FASTA, all annotation
coordinates were in range, and recorded a fresh FAI checksum.

## Executed routes and observations

ONT executed NanoPlot, Pychopper with `PCS109`, minimap2 `-x splice`, streamed
samtools sort/index, shared samtools QC, and normalization. NanoPlot observed
4,000 reads and 1,386,059 bases. Pychopper reported 629 primer-classified reads,
8 rescued reads, and 3,367 unusable reads; 637 classified-plus-rescued reads
entered alignment. One read aligned to the deliberately small reference with a
splice-aware CIGAR containing `N` operations. The others remained in the BAM as
unmapped records.

PacBio HiFi executed NanoPlot `--ubam`, lima `--isoseq --peek-guess`, explicit
selection of `IsoSeqX_bc06_5p--IsoSeqX_3p`, refine with `--require-polya`,
cluster2 with singletons disabled, pbmm2 `--preset ISOSEQ --unmapped`, shared
samtools QC, and normalization. Lima assigned 751 reads to the selected pair;
refine emitted 751 FLNC reads; cluster2 emitted 33 non-singleton consensus
records; pbmm2 retained all 33 records, of which 4 mapped to chrM and 29 were
unmapped.

The `isoseq_flnc` route was executed separately using that genuine post-refine
751-read FLNC BAM. It bypassed lima/refine, reproduced 33 cluster2 records, and
produced the same 4 mapped plus 29 unmapped counts. No reliable universal BAM
signature was identified for distinguishing every third-party FLNC BAM from a
clustered transcript BAM; truly post-refine FLNC data remains an explicit
semantic input contract.

NanoPlot 1.48.0 completed on the PacBio HiFi uBAM and recovered 10,000 reads,
19,250,732 bases, read N50 2,002, and CCS-derived base-quality summaries. Those
qualities are meaningful as BAM base qualities, but ONT pore/channel metrics
are not interpreted for PacBio. Its bundled browser failed to export static PNG
plots in this Docker environment; NanoStats, interactive HTML plots, and the
HTML report were complete and parseable.

For every final BAM, samtools quickcheck, header inspection, idxstats, flagstat,
stats, and coverage completed. The workflow verified coordinate sort order,
readable BAI, non-empty records, exact FASTA/BAM sequence dictionaries, and the
expected sample/platform read group. Every normalized BAM, BAI, QC directory,
runtime-version file, and reference-validation path in each manifest existed.

## Runtime pins actually observed

| Tool | Observed version | Executed image digest |
|---|---|---|
| NanoPlot | 1.48.0 | `sha256:5e790a28aa8ec2974180415c3208e89a1741a5adf4e5de822f83a73dc64d0218` |
| Pychopper | 2.7.10 | `sha256:c311ea485aa6a88cab67542b8cae473d6670f4f33d48d1d467f9d960ee57ff32` |
| minimap2 | 2.30-r1287 | `sha256:a46a2b8c93894d01d1e32ed4fee123c2fc1e9febbe2a513526f19b66c2338c66` |
| samtools (ONT image) | 1.23.1 | `sha256:a46a2b8c93894d01d1e32ed4fee123c2fc1e9febbe2a513526f19b66c2338c66` |
| samtools (shared QC) | 1.24 | `sha256:a55ddea590e567a91df592300a960aa534cfc1bd16e7623e3938ec21f4f3df15` |
| lima | 26.2.1 | `sha256:34594627f817547948800bb9d071b85a57176cb14989e44d782600e78f77e64b` |
| Iso-Seq CLI | 26.2.0 | `sha256:a2936a2c9249e4a5b471501af5fa09704d7356cc98042755b510087ee9a1638f` |
| pbmm2 | 26.2.0 | `sha256:9ad5f7d0295befb8f46799f44955eae12344ebe5a5d3513495c8d563081404ef` |
| IsoQuant | 4.0.0 | `sha256:ff40b81d8945f6724184fb5be90fb3740cbf9b6fe2680db4cbeed816bab73142` |
| SQANTI3 | 6.0.2 | `sha256:3ebeb658fb2e0e308e73cba50dd3b97c63963d8225f303f298ad8ee76264939d` |

## Block 2 real-execution results

Both normalized real samples passed the Block 1-to-Block 2 reference-identity
contract and ran IsoQuant 4.0.0 in bulk mode. ONT used the `nanopore` data type;
PacBio used `pacbio_ccs`; both were marked full-length from their preprocessing
history. The optional IsoQuant canonical-junction annotation was disabled after
the pinned runtime crashed in that reporting path; model construction remained
tool-native.

The compact ONT fixture produced no transcript models under IsoQuant's native
support policy. The workflow retained a comments-only raw GTF, an empty
transcript FASTA, schema-bearing empty classification/junction/structural
tables, a zero-count structural summary, and `NO_TRANSCRIPT_MODELS` SQANTI3
status. No threshold was lowered and no model was fabricated.

The compact PacBio fixture produced one supported transcript model:
`ENST00000361789.2` on `chrM:14747-15887` (`+` strand), 1,141 nt, one exon,
and one supporting read. SQANTI3 classified it as `full-splice_match` against
`ENSG00000198727.2` / `ENST00000361789.2`. The corrected FASTA contained the
same transcript ID and exactly 1,141 sequence bases. The raw GTF,
classification, junction schema, structural-identity TSV, JSON/Markdown
summaries, and contract paths all passed the independent integration checker.

The pinned SQANTI3 image was executed directly before the workflow. Its
observed CLI reported `SQANTI3 6.0.2` and accepted GTF isoforms through
`--isoforms`, the reference through `--refGTF`/`--refFasta`, and the supported
short options `-t`, `-n`, `-o`, and `-d`. IsoQuant-to-SQANTI3 interoperability
was therefore exercised with a real non-empty PacBio model as well as the real
empty-model ONT path.

Per-sample transcriptome contracts use schema version 2.0 and link the raw GTF,
conditional FASTA, native SQANTI3 classification/junction files, IsoQuant
support/output directory, SQANTI3 output directory, structural identity and
summary artifacts, Block 1 manifest, and exact Block 2 version records. The
observed version files record IsoQuant 4.0.0 and SQANTI3 6.0.2 together with
their immutable image digests and execution policies.

## Hardening resulting from execution

- Multiplexed lima output requires `primer_pair` when multiple BAMs appear.
- pbmm2 uses BAM-compatible `--sample`, not the FASTA/FASTQ-only `--rg`.
- pbmm2 includes `--unmapped`, making QC denominators consistent with ONT.
- Reference compatibility and checksum validation precede scientific branches.
- Alignment validation checks BAM, BAI, sorting, dictionaries, read groups, and
  non-empty content.
- Manifests reference aggregated, observed runtime-version artifacts.
- Real-fixture FASTA generation now uses uniformly wrapped LF-only records;
  sequence identity is checked across normalization.
- IsoQuant 4.0.0's crashing optional `--check_canonical` reporting path is not
  enabled; SQANTI3 remains the independent structural classifier.
- SQANTI3 invocation uses only its observed 6.0.2 CLI. A task-local
  `libbz2.so.1` compatibility symlink supplies the loader name omitted by the
  pinned image without modifying the image.
- SQANTI3's optional R HTML report is skipped because version 6.0.2 aborts on a
  valid zero-junction classification. Required native tables, corrected
  GTF/FASTA, and standardized summaries are unaffected.
- Integration checks require exact FASTA/classification ID and length agreement
  and exact parseable Block 2 runtime versions.

Known limitations are the intentionally tiny two-contig reference, different
biological sources for the platform fixtures, environment-sensitive NanoPlot
static PNG export, the skipped SQANTI3 R HTML report, and the inability to infer
biological performance from execution fixtures. The ONT fixture's zero-model
outcome and the PacBio fixture's single model are execution-contract evidence,
not estimates of caller sensitivity or biological composition.

## Block 3 real-execution validation

Block 3 benchmark mode completed on 2026-10-07/08 with the same normalized ONT
and PacBio BAM/BAI records and the same mini-reference/annotation. IsoQuant
reused the Block 2 module. FLAIR and Bambu consumed the normalized genome BAMs
directly; neither received a separate alignment or extra junction evidence.
The first benchmark invocation missed old upstream cache entries and reproduced
Block 1 successfully; later fix-validation runs used Nextflow resume (final
run: 21 cached tasks, 16 successful tasks).

| Tool | Observed version | Container identity |
|---|---:|---|
| FLAIR | 3.0.1 | `quay.io/biocontainers/flair:3.0.1--pyhdfd78af_0@sha256:a39894c1507c68241fba62ddbe1a931d5a6bfaacbcb6be24631d46b60e67ce84` |
| Bambu | 3.14.0 | `longread-rnaseq-pipeline/bambu:3.14.0-bioc3.23-r4.6`, observed local image ID `sha256:8e073d65a08cd3924b4e43491e9b20f191408e1e2c3a4cf4ee064fe0eb33b647` |
| Bioconductor | 3.23 | base image index `sha256:821dbf9ac119eac41f177531c7ca8fc7084c99c04eb218a1aa7f52cb63bad5d9` |
| R | 4.6.1 | same Bambu image |

FLAIR used `transcriptome --genomealignedbam`, its native minimum support of
three, 100 bp caller-native end clustering, and `--noaligntoannot`. Both compact
BAMs lacked enough usable junction evidence; FLAIR inspected each BAM and
reported no correctable junctions. The workflow records
`NO_USABLE_SPLICE_JUNCTIONS` and publishes valid empty GTF/BED/FASTA/read-map
and support artifacts. No support threshold was lowered.

Bambu ran discovery and quantification as separate documented phases. The ONT
fixture produced no read-supported models and triggered Bambu 3.14.0's sparse
empty-equivalence-class quantification error; the valid discovery object is
retained, quantification is explicitly unavailable for sparse input, and empty
supported-model GTF/count artifacts are published. PacBio quantification
completed, but Bambu produced no models satisfying its native `readCount >= 1`
comparison-universe rule. Complete discovery objects and all-model GTFs remain
available; unobserved supplied-reference transcripts are not counted as
reconstructions.

The common SQANTI3 6.0.2 module ran or emitted its explicit empty-model contract
for all six sample/caller combinations. ONT was empty for all three callers.
PacBio retained the one IsoQuant single-exon mitochondrial FSM from Block 2;
FLAIR and Bambu were empty. The observed PacBio membership was therefore one
`isoquant_only` FSM. Multi-exon, intron-chain, and junction unions were empty,
so their Jaccards are `null` with `both_empty`; single-exon comparisons against
the empty challengers have intersection 0, union 1, Jaccard 0, and `one_empty`.
These results validate sparse contracts and do not estimate caller performance.

Execution hardening was limited to new benchmark code: Bambu's sparse
quantification and zero-row GTF serialization are explicit; FLAIR's precise
no-usable-junction outcome maps to an empty caller while other errors remain
fatal; and fixed caller artifact names avoid shell-glob ambiguity. The
independent checker validated caller and SQANTI3 directories, standardized
tables, FASTA/classification IDs, memberships, metrics, summaries, contracts,
published paths, and aggregated runtime provenance. Deterministic fixtures
separately exercise meaningful concordance logic. No consensus transcriptome or
downstream biological analysis was created.

## Block 4 real-execution validation

Block 4 biology mode completed on 2026-10-08 using the production IsoQuant
transcriptome and the same genuine compact ONT and PacBio fixtures. FLAIR and
Bambu benchmark outputs were not used to define or filter the production
models. The final resumed run completed 36 processes (31 cached and 5
executed), and the independent biology-contract checker passed.

The two fixtures retain the common descriptive condition label but originate
from different biological sources and platforms. No contrast was supplied.
Gene DE, transcript DE, DTU, differential splicing, and isoform switching all
therefore emitted explicit `NOT_RUN_NO_EXPLICIT_CONTRAST` statuses; ONT was not
compared biologically with PacBio. Expression QC, abundance matrices,
strand-aware transcript boundaries, integrated evidence, and provenance
completed. The single real model generated no SUPPA event, giving
`COMPLETE_EMPTY_RESULT` for events and `NOT_RUN_NO_EVENTS` for PSI and
differential splicing.

The PacBio path contained one transcript. `TransDecoder.LongOrfs` completed and
retained one 187-aa candidate, while `TransDecoder.Predict` encountered its
tool-native small-training-set/PWM failure. The workflow recorded
`PREDICTION_NOT_AVAILABLE_SPARSE_INPUT`, retained the candidate as
`ORF_CANDIDATE_AVAILABLE`, reported zero predicted-primary transcripts and
`primary_orfs_fabricated: false`, and continued successfully. The transcript
was not labelled non-coding. Empty input remains distinct from
`NO_ORF_CANDIDATE`; unrecognized prediction failures and genuine `LongOrfs`
failures remain fatal.

Observed Block 4 versions were R 4.6.1, Bioconductor 3.23, edgeR 4.10.5,
satuRn 1.20.0, IsoformSwitchAnalyzeR 2.12.0, clusterProfiler 4.20.2,
org.Hs.eg.db 3.23.1, SUPPA release 2.4 (its CLI reports 2.3), and TransDecoder
6.0.0. The contract suite also exercised deterministic replicated statistical
fixtures for edgeR, satuRn, SUPPA, and IsoformSwitchAnalyzeR; those fixtures are
mathematical plumbing tests, not biological validation data. No small real
replicated condition dataset was added, so real biological differential
validation remains pending.

## Block 5 final-report validation

Block 5 reporting completed on 2026-10-08 for both existing compact real-data
routes. The biology report consumed the frozen Block 1, Block 2, and Block 4
contracts; a final resume check reused all 37 scientific/QC tasks and executed
only `FINAL_REPORT`. The benchmark report consumed the Block 1, Block 2, and
Block 3 contracts; its final report-only correction reused all 38 upstream
tasks. The production integrated-biology table retained SHA-256
`3f5dabde2453f03da7a9ef011621ee751b7304bdd4be0a4c4a7b6f3b33e57db6`.

Both modes produced an offline custom HTML report, JSON run summary, TSV
exports, MultiQC report, consolidated provenance, reproducibility manifest,
resolved parameters, and tool/container audit. Independent link checks passed
for every relative report link. The biology report showed the ONT zero-model
case, explicit non-run statistical branches, empty SUPPA event result, and the
PacBio `PREDICTION_NOT_AVAILABLE_SPARSE_INPUT` ORF state without presenting any
of them as zero biological effect. The benchmark report preserved all-empty
set statuses and null Jaccards where unions were empty.

MultiQC 1.35 ran from
`ghcr.io/multiqc/multiqc:v1.35@sha256:976be2094a3bc1dab7315ad03440d43aab09bf7756f84f3a92ef7ae1013e4bfd`.
Its network/version checks and automated summary features were disabled. The
distributed MultiQC HTML bundle contains inactive generic interface code, but
the generated summary wrapper is hidden and no summary request is made.

The deterministic reporting regression suite passed complete, zero-model,
all-empty benchmark, sparse ORF, missing contrast, insufficient replication,
disabled analysis, missing optional output, and failed-required-stage cases.
The CI-equivalent Nextflow stub suite passed, including expected early failures
for invalid configuration. Docker was the validated execution backend.
Apptainer and Singularity profiles resolved through Nextflow configuration, but
neither runtime was installed, so container execution with those backends
remains unvalidated.
