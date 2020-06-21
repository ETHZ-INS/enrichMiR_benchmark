#!/bin/bash

a=/reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes.gtf
fc_opts="-O --largestOverlap --nonOverlap 3 --fracOverlap 0.9 --primary -s 2 -T 2 --extraAttributes gene_name,gene_biotype"
featureCounts -t exon $fc_opts -a $a -o spliced.counts /mnt/schratt/p1006/Cherone_2019/STAR_2/*.sam
featureCounts -t exon --nonSplitOnly $fc_opts -a $a -o exonicNonSplit.counts /mnt/schratt/p1006/Cherone_2019/STAR_2/*.sam
featureCounts -t transcript --nonSplitOnly $fc_opts -a $a -o unspliced.tmp.counts /mnt/schratt/p1006/Cherone_2019/STAR_2/*.sam

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
      colnames(x) <- sapply(strsplit(colnames(x),"\\\."), function(y) rev(y)[1])
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

sp <- fc2se("./spliced.counts")
un <- fc2se("./unspliced.tmp.counts")
tmp <- fc2se("./exonicNonSplit.counts")
rowData(sp)[["Length.unspliced"]] <- rowData(un)[["Length"]]
assayNames(sp) <- "spliced"
tmp <- assay(un)-assay(tmp)
tmp[tmp<0] <- 0
assays(sp)[["unspliced"]] <- tmp
saveRDS(sp, file="cherone.splicing.raw.SE.rds")

RSCRIPT

