#!/usr/bin/env Rscript
suppressPackageStartupMessages(library(IsoformSwitchAnalyzeR))
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=6) stop("usage: run_isoform_switch.R DATASET OUTPUT NUMERATOR DENOMINATOR FDR DIF")
dataset<-args[[1]]; output<-args[[2]]; numerator<-args[[3]]; denominator<-args[[4]]; alpha<-as.numeric(args[[5]]); dif<-as.numeric(args[[6]])
dir.create(output,recursive=TRUE,showWarnings=FALSE)
design_status<-jsonlite::fromJSON(file.path(dataset,"design_status.json")); meta<-read.delim(file.path(dataset,"sample_metadata.tsv"),stringsAsFactors=FALSE)
empty<-data.frame(gene_id=character(),transcript_id=character(),dIF=numeric(),p_value=numeric(),adjusted_p_value=numeric(),condition_1=character(),condition_2=character(),statistical_switch=logical(),functional_consequence_evidence=character())
write.table(empty,file.path(output,"isoform_switches.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(),file.path(output,"predicted_consequences.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write_status<-function(s,details=list()) jsonlite::write_json(c(list(status=s,engine="IsoformSwitchAnalyzeR",test_engine="satuRn",alpha=alpha,dIF_cutoff=dif),details),file.path(output,"isoform_switch_status.json"),pretty=TRUE,auto_unbox=TRUE)
if(design_status$status!="READY") { write_status(design_status$status,list(reasons=design_status$reasons))
} else {
  counts<-read.delim(file.path(dataset,"transcript_counts.tsv"),check.names=FALSE); names(counts)[1]<-"isoform_id"
  expression<-read.delim(file.path(dataset,"transcript_tpm.tsv"),check.names=FALSE); names(expression)[1]<-"isoform_id"
  design<-data.frame(sampleID=meta$sample,condition=meta$condition,stringsAsFactors=FALSE)
  comparisons<-data.frame(condition_1=denominator,condition_2=numerator)
  result<-tryCatch({
    switch<-importRdata(isoformCountMatrix=counts,isoformRepExpression=expression,designMatrix=design,isoformExonAnnoation=file.path(dataset,"production_transcripts.gtf"),isoformNtFasta=file.path(dataset,"production_transcripts.fasta"),comparisonsToMake=comparisons,detectUnwantedEffects=FALSE,addAnnotatedORFs=FALSE,fixStringTieAnnotationProblem=FALSE,estimateDifferentialGeneRange=FALSE,showProgress=FALSE,quiet=TRUE)
    switch<-preFilter(switch,isoCount=0,min.Count.prop=0,IFcutoff=0,min.IF.prop=0,removeSingleIsoformGenes=TRUE,quiet=TRUE)
    switch<-isoformSwitchTestSatuRn(switch,alpha=alpha,dIFcutoff=dif,reduceToSwitchingGenes=FALSE,reduceFurtherToGenesWithConsequencePotential=FALSE,diagplots=FALSE,showProgress=FALSE,quiet=TRUE)
    native<-switch$isoformSwitchAnalysis
    writeLines(names(native),file.path(output,"native_result_columns.txt"))
    capture.output(str(native),file=file.path(output,"native_result_structure.txt"))
    features<-switch$isoformFeatures
    writeLines(names(features),file.path(output,"native_feature_columns.txt"))
    capture.output(str(features),file=file.path(output,"native_feature_structure.txt"))
    native<-merge(native,features[,intersect(c("iso_ref","gene_id","dIF"),names(features)),drop=FALSE],by="iso_ref",all.x=TRUE,sort=FALSE)
    out<-data.frame(gene_id=native$gene_id,transcript_id=native$isoform_id,dIF=native$dIF,p_value=native$pval,adjusted_p_value=native$padj,condition_1=native$condition_1,condition_2=native$condition_2,statistical_switch=native$padj<alpha & abs(native$dIF)>=dif,functional_consequence_evidence="NOT_EVALUATED_WITHOUT_OPTIONAL_DOMAIN_HOMOLOGY_DATA")
    write.table(out,file.path(output,"isoform_switches.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    saveRDS(switch,file.path(output,"isoform_switch_analysis.rds")); write_status(if(nrow(out)) "COMPLETE" else "COMPLETE_EMPTY_RESULT",list(transcripts_tested=nrow(out))); TRUE
  },error=function(e){write_status("NOT_RUN_TOOL_ERROR",list(error=conditionMessage(e))); FALSE})
}
jsonlite::write_json(list(R=as.character(getRversion()),Bioconductor=as.character(BiocManager::version()),IsoformSwitchAnalyzeR=as.character(packageVersion("IsoformSwitchAnalyzeR"))),file.path(output,"versions.json"),pretty=TRUE,auto_unbox=TRUE)
