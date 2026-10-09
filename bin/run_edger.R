#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(edgeR))
args <- commandArgs(trailingOnly=TRUE)
if (length(args) != 6) stop("usage: run_edger.R DATASET OUTPUT NUMERATOR DENOMINATOR FDR DESIGN_FORMULA")
dataset <- args[[1]]; output <- args[[2]]; numerator <- args[[3]]; denominator <- args[[4]]; fdr_cutoff <- as.numeric(args[[5]]); design_formula <- args[[6]]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
design_status <- jsonlite::fromJSON(file.path(dataset, "design_status.json"))
metadata <- read.delim(file.path(dataset, "sample_metadata.tsv"), check.names=FALSE, stringsAsFactors=FALSE)

empty_result <- function(path, feature_name) {
  write.table(setNames(data.frame(character(), numeric(), numeric(), numeric(), numeric(), numeric(), character()),
    c(feature_name,"logFC","logCPM","F","PValue","FDR","contrast")), path, sep="\t", quote=FALSE, row.names=FALSE)
}
write_status <- function(path, status, details=list()) jsonlite::write_json(c(list(status=status, engine="edgeR", workflow="quasi-likelihood GLM", contrast=paste0(numerator,"-",denominator), fdr_cutoff=fdr_cutoff), details), path, pretty=TRUE, auto_unbox=TRUE)

run_level <- function(level) {
  feature_name <- if (level == "genes") "gene_id" else "transcript_id"
  result_path <- file.path(output, paste0(level, ".differential_expression.tsv"))
  filter_path <- file.path(output, paste0(level, ".filtering.tsv"))
  status_path <- file.path(output, paste0(level, ".status.json"))
  empty_result(result_path, feature_name)
  x <- read.delim(file.path(dataset, if(level == "genes") "gene_counts.tsv" else "transcript_counts.tsv"), check.names=FALSE)
  if (design_status$status != "READY") { write.table(data.frame(feature_id=x[[1]], tested=FALSE, reason=design_status$status), filter_path, sep="\t", quote=FALSE, row.names=FALSE); write_status(status_path, design_status$status, list(reasons=design_status$reasons)); return() }
  counts <- as.matrix(x[,-1,drop=FALSE]); rownames(counts) <- x[[1]]; storage.mode(counts) <- "integer"
  if (!nrow(counts) || all(rowSums(counts)==0)) { write_status(status_path, "NOT_RUN_ALL_ZERO_FEATURES"); return() }
  metadata$condition <- factor(metadata$condition)
  design <- tryCatch(model.matrix(as.formula(design_formula), metadata), error=function(e) NULL)
  if (is.null(design)) { write_status(status_path, "NOT_RUN_INVALID_DESIGN_FORMULA"); return() }
  colnames(design) <- sub("^condition", "", colnames(design))
  if (qr(design)$rank < ncol(design)) { write_status(status_path, "NOT_RUN_CONFOUNDED_DESIGN"); return() }
  y <- DGEList(counts=counts); keep <- filterByExpr(y, design=design)
  write.table(data.frame(feature_id=rownames(y), tested=keep, reason=ifelse(keep,"TESTED","INSUFFICIENT_EXPRESSION")), filter_path, sep="\t", quote=FALSE, row.names=FALSE)
  if (!any(keep)) { write_status(status_path, "NOT_RUN_INSUFFICIENT_EXPRESSION"); return() }
  y <- y[keep,,keep.lib.sizes=FALSE]; y <- calcNormFactors(y); y <- estimateDisp(y, design, robust=TRUE)
  fit <- glmQLFit(y, design, robust=TRUE)
  if (!(numerator %in% colnames(design)) || !(denominator %in% colnames(design))) { write_status(status_path, "NOT_RUN_INVALID_CONTRAST"); return() }
  contrast <- rep(0, ncol(design)); names(contrast) <- colnames(design); contrast[numerator] <- 1; contrast[denominator] <- -1
  test <- glmQLFTest(fit, contrast=contrast); result <- topTags(test, n=Inf, sort.by="none")$table
  result[[feature_name]] <- rownames(result); result$contrast <- paste0(numerator,"-",denominator)
  result <- result[,c(feature_name,"logFC","logCPM","F","PValue","FDR","contrast")]
  if (level == "transcripts") {
    tm <- read.delim(file.path(dataset,"transcript_metadata.tsv"), stringsAsFactors=FALSE)
    result <- merge(result, tm[,intersect(c("transcript_id","gene_id","structural_category","known_novel"),names(tm)),drop=FALSE], by="transcript_id", all.x=TRUE, sort=FALSE)
  }
  write.table(result, result_path, sep="\t", quote=FALSE, row.names=FALSE, na="NA")
  write_status(status_path, if(nrow(result)) "COMPLETE" else "COMPLETE_EMPTY_RESULT", list(features_tested=nrow(result), significant_fdr=sum(result$FDR <= fdr_cutoff, na.rm=TRUE), filter="edgeR::filterByExpr"))
}
run_level("genes"); run_level("transcripts")
versions <- list(R=as.character(getRversion()), Bioconductor=as.character(BiocManager::version()), edgeR=as.character(packageVersion("edgeR")))
jsonlite::write_json(versions, file.path(output,"versions.json"), pretty=TRUE, auto_unbox=TRUE)
