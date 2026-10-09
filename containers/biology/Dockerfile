FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:821dbf9ac119eac41f177531c7ca8fc7084c99c04eb218a1aa7f52cb63bad5d9

RUN Rscript -e 'options(repos = BiocManager::repositories()); BiocManager::install(c("edgeR", "satuRn", "IsoformSwitchAnalyzeR", "clusterProfiler", "org.Hs.eg.db", "jsonlite"), ask = FALSE, update = FALSE, Ncpus = 4)'
RUN Rscript -e 'stopifnot(as.character(packageVersion("edgeR")) == "4.10.5", as.character(packageVersion("satuRn")) == "1.20.0", as.character(packageVersion("IsoformSwitchAnalyzeR")) == "2.12.0", as.character(packageVersion("clusterProfiler")) == "4.20.2", as.character(packageVersion("org.Hs.eg.db")) == "3.23.1")'

LABEL org.opencontainers.image.title="longread-rnaseq-pipeline Block 4 statistics" \
      org.opencontainers.image.version="bioc3.23-r4.6"
