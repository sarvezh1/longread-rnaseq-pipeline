#!/usr/bin/env Rscript
suppressPackageStartupMessages({library(clusterProfiler); library(org.Hs.eg.db); library(AnnotationDbi)})
args<-commandArgs(trailingOnly=TRUE)
if(length(args)!=5) stop("usage: run_enrichment.R DATASET STATISTICS DTU SWITCHES OUTPUT")
dataset<-args[[1]]; statistics<-args[[2]]; dtu_dir<-args[[3]]; switch_dir<-args[[4]]; output<-args[[5]]
dir.create(output,recursive=TRUE,showWarnings=FALSE)
enabled<-tolower(Sys.getenv("BLOCK4_ENABLE_ENRICHMENT","false"))=="true"; fdr<-as.numeric(Sys.getenv("BLOCK4_FDR","0.05"))
write.table(data.frame(gene_set=character(),ontology=character(),ID=character(),Description=character(),GeneRatio=character(),BgRatio=character(),pvalue=numeric(),p.adjust=numeric(),qvalue=numeric(),geneID=character(),Count=integer()),file.path(output,"go_enrichment.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
write.table(data.frame(input_gene_id=character(),normalized_ensembl_id=character(),entrez_id=character(),mapped=logical()),file.path(output,"identifier_mapping.tsv"),sep="\t",quote=FALSE,row.names=FALSE)
status<-function(s,d=list()) jsonlite::write_json(c(list(status=s,engine="clusterProfiler",database="org.Hs.eg.db",background_policy="testable feature universe",automated_natural_language_interpretation=FALSE),d),file.path(output,"enrichment_status.json"),pretty=TRUE,auto_unbox=TRUE)
if(!enabled){status("NOT_RUN_DISABLED")
} else {
  de<-read.delim(file.path(statistics,"genes.differential_expression.tsv"),stringsAsFactors=FALSE)
  universe<-unique(sub("\\..*$","",de$gene_id)); selected<-unique(sub("\\..*$","",de$gene_id[!is.na(de$FDR)&de$FDR<=fdr]))
  if(!length(selected)){status("NOT_RUN_EMPTY_GENE_SET",list(testable_background_size=length(universe)))
  } else {
    mapping<-AnnotationDbi::select(org.Hs.eg.db,keys=unique(c(universe,selected)),keytype="ENSEMBL",columns="ENTREZID")
    mapout<-data.frame(input_gene_id=mapping$ENSEMBL,normalized_ensembl_id=mapping$ENSEMBL,entrez_id=mapping$ENTREZID,mapped=!is.na(mapping$ENTREZID)); write.table(mapout,file.path(output,"identifier_mapping.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")
    sel<-unique(mapping$ENTREZID[mapping$ENSEMBL%in%selected & !is.na(mapping$ENTREZID)]); bg<-unique(mapping$ENTREZID[mapping$ENSEMBL%in%universe & !is.na(mapping$ENTREZID)])
    if(!length(sel)){status("NOT_RUN_NO_MAPPED_SELECTED_GENES",list(input_selected=length(selected),unmapped_selected=length(selected)))
    } else {
      result<-enrichGO(gene=sel,universe=bg,OrgDb=org.Hs.eg.db,keyType="ENTREZID",ont="ALL",pAdjustMethod="BH",readable=FALSE)
      out<-as.data.frame(result); if(nrow(out)){out$gene_set<-"differentially_expressed_genes"; out$ontology<-out$ONTOLOGY; out<-out[,c("gene_set","ontology","ID","Description","GeneRatio","BgRatio","pvalue","p.adjust","qvalue","geneID","Count")]; write.table(out,file.path(output,"go_enrichment.tsv"),sep="\t",quote=FALSE,row.names=FALSE,na="NA")}
      status(if(nrow(out))"COMPLETE" else "COMPLETE_EMPTY_RESULT",list(input_selected=length(selected),mapped_selected=length(sel),testable_background_size=length(universe),mapped_background_size=length(bg)))
    }
  }
}
jsonlite::write_json(list(R=as.character(getRversion()),Bioconductor=as.character(BiocManager::version()),clusterProfiler=as.character(packageVersion("clusterProfiler")),org.Hs.eg.db=as.character(packageVersion("org.Hs.eg.db"))),file.path(output,"versions.json"),pretty=TRUE,auto_unbox=TRUE)
