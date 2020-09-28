devtools::load_all("/mnt/schratt/tgermade_test/master_19_20/enrichMiR/enrichMiR/")

cleanDEA <- function(dea){
  
  # find correct column names & order
  fdr.n <- if("padj" %in% colnames(dea)) "padj" else if("FDR" %in% colnames(dea)) "FDR"
  pval.n <- if("pvalue" %in% colnames(dea)) "pvalue" else if("PValue" %in% colnames(dea)) "PValue"
  dea <- dea[order(dea[,fdr.n], dea[,pval.n]),]
  # when using gene-level annotations: clean up rownames & aggregate tx
  if(length(grep("^ENSG|^ENSRNOG|^ENSMUSG", row.names(dea))) > 100){
    g <- sapply(strsplit(row.names(dea), "\\."), FUN=function(x) x[2])
    dea <- dea[!duplicated(g),]
    row.names(dea) <- g[!duplicated(g)]
  }
  # when using transcript-level annotations: remove transcript version
  if(length(grep("^ENST", rownames(dea))) > 100){
    g <- sapply(strsplit(row.names(dea), "\\."), FUN=function(x) x[1])
    row.names(dea) <- g[!duplicated(g)]
  }
  # get rid of Infinite values
  for(i in 1:ncol(dea)){
    if(any(is.infinite(dea[,i]))){
      w <- which(is.infinite(dea[,i]))
      dea[w,i] <- max(abs(dea[setdiff(which(dea[,fdr.n]<0.5),w),i]))*sign(dea[w,i])
    }
  }
  # get rid of NAs
  dea[is.na(dea[,pval.n]),pval.n] <- 1
  dea[is.na(dea[,fdr.n]),fdr.n] <- 1
  
  return(dea)
}

suppressPackageStartupMessages(c(
  library(SummarizedExperiment)
))


datasets <- c("amin","jeong")
for(ds in datasets){

  if(ds=="amin"){
    DEA <- readRDS("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/amin.DEA.SE.rds")
    TS <- readRDS("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/TargetScan_mouse.rds")
  } else {
    DEA <- readRDS("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/jeong.DEA.SE.rds")
    TS <- readRDS("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/TargetScan_human.rds")
  }
  DEA <- rowData(DEA)$DEA.all
  DEA <- cleanDEA(DEA)
  if(all(c("family", "feature") %in% colnames(TS))){
    colnames(TS) <- gsub("^family$","set", colnames(TS))
  }
  TS <- DataFrame(TS)
  
  # measure time for each test (in seconds)
  tests <- c("overlap","siteoverlap","woverlap","modscore","modsites","regmir","regmirb","ks","mw","gsea","areamir")
  tests <- setdiff(tests,"gsea")
  t <- NULL
  e <- NULL
  for(test in tests){
    start_time <- Sys.time()
    e <- list(e, enrichMiR(DEA=as.data.frame(DEA), TS=TS, miRNA.expression=NULL, families=metadata(TS)$families, cleanNames=TRUE, tests=test) )
    end_time <- Sys.time()
    t <- c(t, end_time - start_time)
  }
  names(t) <- tests
  e <- unlist(e)
  names(e) <- tests
  
  if(ds=="amin"){
    saveRDS(t, "results/amin.runtimes.rds")
  } else {
    saveRDS(t, "results/jeong.runtimes.rds")
  }
}



