#!/usr/bin/env Rscript

suppressPackageStartupMessages(library(bambu))
suppressPackageStartupMessages(library(SummarizedExperiment))

arguments <- commandArgs(trailingOnly = TRUE)
if (length(arguments) != 5) {
    stop("usage: run_bambu.R BAM BAI GTF FASTA OUTPUT_DIR")
}
bam <- normalizePath(arguments[[1]])
bai <- normalizePath(arguments[[2]])
gtf <- normalizePath(arguments[[3]])
fasta <- normalizePath(arguments[[4]])
output_dir <- arguments[[5]]
dir.create(output_dir, recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "all_models"), recursive = TRUE, showWarnings = FALSE)
dir.create(file.path(output_dir, "detected_models"), recursive = TRUE, showWarnings = FALSE)

write_gtf_or_empty <- function(models, path, label) {
    if (length(models) == 0) {
        writeLines(paste0("# no ", label, " transcript models"), path)
    } else {
        writeToGTF(models, path)
    }
}

# The .bai is an explicit workflow input even though Rsamtools discovers it by BAM basename.
if (!file.exists(bai)) stop("BAM index is missing")
annotation <- prepareAnnotations(gtf)
discovered <- bambu(
    reads = bam,
    annotations = annotation,
    genome = fasta,
    discovery = TRUE,
    quant = FALSE,
    ncore = as.integer(Sys.getenv("NXF_TASK_CPUS", "1"))
)
saveRDS(discovered, file.path(output_dir, "bambu_discovery.rds"))
write_gtf_or_empty(discovered, file.path(output_dir, "all_models", "all_models_extended_annotations.gtf"), "Bambu")

# The discovery object includes the supplied annotation. Compare only models for which Bambu
# reports native read support; the complete object is retained above. This is not a relaxed
# discovery threshold. It prevents unobserved reference transcripts entering the benchmark.
read_count <- if ("readCount" %in% colnames(mcols(discovered))) mcols(discovered)$readCount else rep(NA_real_, length(discovered))
keep <- !is.na(read_count) & read_count >= 1
detected <- discovered[keep]
write_gtf_or_empty(detected, file.path(output_dir, "detected_models", "detected_extended_annotations.gtf"), "supported Bambu")

# Quantification is attempted against the discovered annotation as documented by Bambu. Very
# sparse inputs can trigger an empty-equivalence-class type error in Bambu 3.14.0; discovery
# remains valid and the missing abundance is made explicit rather than weakening thresholds.
quantification_error <- NULL
se <- tryCatch(
    bambu(
        reads = bam, annotations = discovered, genome = fasta,
        discovery = FALSE, quant = TRUE,
        ncore = as.integer(Sys.getenv("NXF_TASK_CPUS", "1"))
    ),
    error = function(error) {
        quantification_error <<- conditionMessage(error)
        NULL
    }
)
count_path <- file.path(output_dir, "detected_models", "bambu_transcript_counts.tsv")
if (is.null(se)) {
    write.table(data.frame(transcript_id = names(detected), count = rep(NA_real_, length(detected))), count_path,
                sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines(c("status\tnot_available_sparse_input", paste0("error\t", quantification_error)),
               file.path(output_dir, "quantification_status.tsv"))
} else {
    saveRDS(se, file.path(output_dir, "bambu_result.rds"))
    quantification_dir <- file.path(output_dir, "quantification")
    dir.create(quantification_dir, recursive = TRUE, showWarnings = FALSE)
    writeBambuOutput(se, path = quantification_dir, prefix = "quantified_")
    counts <- assays(se)$counts
    count_table <- data.frame(transcript_id = rownames(counts), count = rowSums(counts, na.rm = TRUE))
    write.table(count_table[count_table$transcript_id %in% names(detected), ], count_path,
                sep = "\t", quote = FALSE, row.names = FALSE)
    writeLines("status\tcomplete", file.path(output_dir, "quantification_status.tsv"))
}

metadata <- data.frame(
    transcript_id = names(discovered),
    novel_transcript = if ("novelTranscript" %in% colnames(mcols(discovered))) mcols(discovered)$novelTranscript else NA,
    read_count = read_count,
    detected_for_comparison = keep,
    stringsAsFactors = FALSE
)
write.table(metadata, file.path(output_dir, "caller_native_classification.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
