#!/bin/bash

# gene annotation file (GTF)
a=/reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_98-2019-09-3/genes.gtf
# featureCounts options
fc_opts="-O --largestOverlap --nonOverlap 3 --fracOverlap 0.9 --primary -s 2 -T 14 --extraAttributes gene_name,gene_biotype"
# create featureCounts files
featureCounts -t exon $fc_opts -a $a -o spliced.counts ../STAR/*Aligned.out.sam
featureCounts -t exon --nonSplitOnly $fc_opts -a $a -o exonicNonSplit.counts ../STAR/*Aligned.out.sam
featureCounts -t transcript --nonSplitOnly $fc_opts -a $a -o unspliced.tmp.counts ../STAR/*Aligned.out.sam

# R script to differentiatie btw. exonic & intronic counts
R --slave --no-restore --no-save <<RSCRIPT

fc2se <- function(countfile){
	library(SummarizedExperiment)
	x <- read.delim(countfile,header=TRUE,skip=1,row.names=1)
	rd <- x[,1:7]
	levels(rd[["Chr"]]) <- sapply(strsplit(levels(rd[["Chr"]]),";"), FUN=function(x){ paste(unique(x), collapse=";")})
	levels(rd[["Strand"]]) <- sapply(strsplit(levels(rd[["Strand"]]),";"), FUN=function(x) x[1])
	x <- as.matrix(x[,-1:-7])
	colnames(x) <- gsub("Aligned.out.sam","",colnames(x),fixed=TRUE)
	colnames(x) <- gsub("...STAR.","",colnames(x),fixed=TRUE)
	se <- SummarizedExperiment(list(counts=x), rowData=rd)
	row.names(se) <- paste(row.names(se),rowData(se)[["gene_name"]],sep=".")
	se[["Assigned"]] <- colSums(x)
	if(file.exists(paste0(countfile,".summary"))){
		x <- read.delim(paste0(countfile,".summary"),header=TRUE,row.names=1)
		colnames(x) <- gsub("Aligned.out.sam","",colnames(x),fixed=TRUE)
		colnames(x) <- gsub("...STAR.","",colnames(x),fixed=TRUE)
		se[["Unassigned"]] <- colSums(x[-1,colnames(se)])
	}
	se
}

sp <- fc2se("./spliced.counts")
un <- fc2se("./unspliced.tmp.counts")
tmp <- fc2se("./exonicNonSplit.counts")
rowData(sp)[["Length.unspliced"]] <- rowData(un)[["Length"]]
assayNames(sp) <- "spliced"
tmp <- assay(un)-assay(tmp)
tmp[tmp<0] <- 0
assays(sp)[["unspliced"]] <- tmp
saveRDS(sp, file="bartel.intron.raw.SE.rds")

RSCRIPT