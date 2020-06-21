#!/bin/bash/

a=/reference/Rattus_norvegicus/Ensembl/Rnor_6.0/Annotation/Genes/genes.gtf
fc_opts=" --largestOverlap --primary -s 2 -T 10 --extraAttributes gene_name,gene_biotype"
featureCounts -t transcript $fc_opts -a $a -o comb.counts /mnt/schratt/p3153_ratPolyA_mirDuplexes/STAR_40416_2019-10-24--09-03-17/*.bam

R --slave --no-restore --no-save <<RSCRIPT

fc2se <- function(countfile){
	library(SummarizedExperiment)
	x <- read.delim(countfile,header=TRUE,skip=1,row.names=1)
	rd <- x[,1:7]
	levels(rd[["Chr"]]) <- sapply(strsplit(levels(rd[["Chr"]]),";"), FUN=function(x){ paste(unique(x), collapse=";")})
	levels(rd[["Strand"]]) <- sapply(strsplit(levels(rd[["Strand"]]),";"), FUN=function(x) x[1])
	x <- as.matrix(x[,-1:-7])
	colnames(x) <- sapply(strsplit(colnames(x),"\\\."), function(y) if(!any(grepl("Neg",y))) paste(rev(y)[[3]],rev(y)[[2]],sep=".") else rev(y)[[2]])
	se <- SummarizedExperiment(list(counts=x), rowData=rd)
	row.names(se) <- paste(row.names(se),rowData(se)[["gene_name"]],sep=".")
	se[["Assigned"]] <- colSums(x)
	if(file.exists(paste0(countfile,".summary"))){
            x <- read.delim(paste0(countfile,".summary"),header=TRUE,row.names=1)
            colnames(x) <- colnames(se)
            se[["Unassigned"]] <- colSums(x[-1,colnames(se)])
	}
	se
}

sp <- fc2se("./comb.counts")
saveRDS(sp, file="ratPolyA.raw.SE.rds")

RSCRIPT
