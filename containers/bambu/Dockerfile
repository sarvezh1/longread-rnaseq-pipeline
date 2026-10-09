FROM bioconductor/bioconductor_docker:RELEASE_3_23@sha256:821dbf9ac119eac41f177531c7ca8fc7084c99c04eb218a1aa7f52cb63bad5d9

ARG BAMBU_VERSION=3.14.0
RUN Rscript -e 'options(repos = BiocManager::repositories()); BiocManager::install("bambu", version = "3.23", ask = FALSE, update = FALSE); stopifnot(as.character(packageVersion("bambu")) == Sys.getenv("BAMBU_VERSION", unset = "3.14.0"))'

LABEL org.opencontainers.image.title="longread-rnaseq-pipeline Bambu runtime" \
      org.opencontainers.image.version="3.14.0-bioc3.23-r4.6" \
      org.opencontainers.image.source="https://bioconductor.org/packages/3.23/bioc/html/bambu.html"
