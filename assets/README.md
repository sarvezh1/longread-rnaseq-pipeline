# Assets

`samplesheet.example.csv` demonstrates every v1 input route. The
`samplesheets/` directory provides focused ONT, PacBio, mixed-platform, and
replicated-condition examples. Their relative paths are illustrative and are
not bundled test data.

Reference genomes, annotations, primers, and sequencing data must not be
committed here. Primer FASTA files are experimental inputs and must be supplied
explicitly by the user; the workflow does not embed a universal primer set.

`reporting/` contains the deterministic status vocabulary, software registry,
MultiQC configuration, and an explicit placeholder used to compose optional
contract channels.
