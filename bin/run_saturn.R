#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(satuRn); library(SummarizedExperiment)})
args <- commandArgs(trailingOnly=TRUE)
if(length(args)!=6) stop("usage: run_saturn.R DATASET OUTPUT NUMERATOR DENOMINATOR FDR DESIGN_FORMULA")
dataset<-args[[1]]; output<-args[[2]]; numerator<-args[[3]]; denominator<-args[[4]]; fdr<-as.numeric(args[[5]]); design_formula<-args[[6]]
dir.create(output, recursive=TRUE, showWarnings=FALSE)
status <- jsonlite::fromJSON(file.path(dataset,"design_status.json")); meta<-read.delim(file.path(dataset,"sample_metadata.tsv"),check.names=FALSE,stringsAsFactors=FALSE)
x<-read.delim(file.path(dataset,"transcript_counts.tsv"),check.names=FALSE); txmeta<-read.delim(file.path(dataset,"transcript_metadata.tsv"),stringsAsFactors=FALSE)
counts<-as.matrix(x[,-1,drop=FALSE]); rownames(counts)<-x[[1]]; storage.mode(counts)<-"integer"
mapping<-txmeta[match(rownames(counts),txmeta$transcript_id),c("transcript_id","gene_id")]
gene_sizes<-table(mapping$gene_id); testable<-!is.na(mapping$gene_id) & gene_sizes[mapping$gene_id]>=2 & rowSums(counts)>0
write.table(data.frame(transcript_id=rownames(counts),gene_id=mapping$gene_id,testable=testable,reason=ifelse(is.na(mapping$gene_id),"MISSING_GENE",ifelse(gene_sizes[mapping$gene_id]<2,"GENE_HAS_ONE_TRANSCRIPT",ifelse(rowSums(counts)==0,"ALL_ZERO","TESTED")))),file.path(output,"dtu_filtering.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
usage<-counts
for(j in seq_len(ncol(counts))) { totals<-rowsum(counts[,j],mapping$gene_id,reorder=FALSE); denom<-totals[match(mapping$gene_id,rownames(totals)),1]; usage[,j]<-ifelse(denom>0,counts[,j]/denom,NA_real_) }
write.table(data.frame(transcript_id=rownames(usage),gene_id=mapping$gene_id,usage,check.names=FALSE),file.path(output,"transcript_usage.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
empty<-data.frame(transcript_id=character(),gene_id=character(),estimate=numeric(),p_value=numeric(),adjusted_p_value=numeric(),contrast=character())
write.table(empty,file.path(output,"dtu_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write_status<-function(s,details=list()) jsonlite::write_json(c(list(status=s,engine="satuRn",contrast=paste0(numerator,"-",denominator),fdr_cutoff=fdr),details),file.path(output,"dtu_status.json"),pretty=TRUE,auto_unbox=TRUE)
if(status$status!="READY") { write_status(status$status,list(reasons=status$reasons))
} else if(sum(testable)<2) { write_status("NOT_RUN_NO_TESTABLE_MULTI_TRANSCRIPT_GENES")
} else {
  se<-SummarizedExperiment(assays=list(counts=counts[testable,,drop=FALSE]),rowData=DataFrame(isoform_id=rownames(counts)[testable],gene_id=mapping$gene_id[testable]),colData=DataFrame(meta,row.names=meta$sample))
  formula<-as.formula(design_formula); design<-model.matrix(formula,meta); colnames(design)<-sub("^condition","",colnames(design))
  if(qr(design)$rank<ncol(design)) { write_status("NOT_RUN_CONFOUNDED_DESIGN")
  } else if(!(numerator %in% colnames(design)) || !(denominator %in% colnames(design))) { write_status("NOT_RUN_INVALID_CONTRAST")
  } else {
  se<-fitDTU(se,formula=formula,parallel=FALSE,verbose=FALSE)
  L<-matrix(0,nrow=ncol(design),ncol=1,dimnames=list(colnames(design),paste0(numerator,"-",denominator))); L[c(numerator,denominator),1]<-c(1,-1)
  se<-testDTU(se,contrasts=L,diagplot1=FALSE,diagplot2=FALSE,sort=FALSE)
  result_name<-grep("fitDTUResult",names(rowData(se)),value=TRUE)[1]; native<-as.data.frame(rowData(se)[[result_name]])
  adjusted <- if("empirical_FDR" %in% names(native) && !all(is.na(native$empirical_FDR))) native$empirical_FDR else native$regular_FDR
  out<-data.frame(transcript_id=rownames(se),gene_id=rowData(se)$gene_id,estimate=native$estimates,p_value=native$pval,adjusted_p_value=adjusted,contrast=paste0(numerator,"-",denominator))
  write.table(out,file.path(output,"dtu_results.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
  write_status(if(nrow(out)) "COMPLETE" else "COMPLETE_EMPTY_RESULT",list(transcripts_tested=nrow(out),significant_fdr=sum(out$adjusted_p_value<=fdr,na.rm=TRUE)))
  }
}
jsonlite::write_json(list(R=as.character(getRversion()),Bioconductor=as.character(BiocManager::version()),satuRn=as.character(packageVersion("satuRn"))),file.path(output,"versions.json"),pretty=TRUE,auto_unbox=TRUE)
