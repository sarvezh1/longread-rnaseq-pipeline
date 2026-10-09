# Site profiles

This directory is reserved for institution- or scheduler-specific Nextflow
configuration. Portable `local`, `docker`, `apptainer`, `singularity`, and
`test` profiles are declared in `nextflow.config`.

The adjacent example configurations show generic scheduler and Apptainer cache
hooks. Copy and adapt one outside the source tree, then supply it with `-c`.
They have not been scheduler-tested. Site-specific queues, accounts, reference
paths, mount options, and credentials must not be committed here.
