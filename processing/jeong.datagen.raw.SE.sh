#!/bin/bash/

a=/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_98-2019-09-3/genes.gtf
fc_opts=" --largestOverlap --primary -s 2 -T 10 --extraAttributes gene_name,gene_biotype"
featureCounts -t transcript $fc_opts -a $a -o comb.counts /mnt/schratt/p1006/Jeong_2017/STAR/*.sam

R --slave --no-restore --no-save <<RSCRIPT

fc2se <- function(countfile){
	library(SummarizedExperiment)
	x <- read.delim(countfile,header=TRUE,skip=1,row.names=1)
	rd <- x[,1:7]
	levels(rd[["Chr"]]) <- sapply(strsplit(levels(rd[["Chr"]]),";"), FUN=function(x){ paste(unique(x), collapse=";")})
	levels(rd[["Strand"]]) <- sapply(strsplit(levels(rd[["Strand"]]),";"), FUN=function(x) x[1])
	x <- as.matrix(x[,-1:-7])
	colnames(x) <- sapply(strsplit(colnames(x),"\\\."), function(y) rev(y)[[3]])
	colnames(x) <- gsub("Aligned","",colnames(x),fixed=TRUE)
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

se <- fc2se("./comb.counts")
saveRDS(se, file="jeong.raw.SE.rds")

RSCRIPT
