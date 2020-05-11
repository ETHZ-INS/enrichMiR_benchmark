suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SEtools)
  library(ggplot2)
  library(plgINS)
  library(edgeR)
  library(DT)
  library(viridis)
  library(RColorBrewer)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(cowplot)
  library(matrixStats)
  library(RUVSeq)
  library(BiocParallel)
})
devtools::load_all("/mnt/schratt/enrichMiR/")

##########################################################################################
# functions

#################################################
#' getDEA
#'
#' @param dea.df dataframe of DEA results (an SE rowData column if generated via DEA() function)
#' @param cleanRowNames set to TRUE if dea.df rownames have ensembl code prefix
#'
#' @return FDR-filtered dataframe containing c("symbol","logFC","PValue","FDR") columns
#'
getDEA <- function(dea.df, cleanRowNames=FALSE){
  if(cleanRowNames){
    symbol <- sapply(rownames(dea.df), function(x) paste(unlist(strsplit(x, split="\\."))[-1], collapse=".") )
    dea.df$symbol <- symbol 
  }
  if(!any(colnames(dea.df) %in% "symbol")) dea.df$symbol <- rownames(dea.df)
  dea <- as.data.frame(dea.df[,c("symbol","logFC", "logCPM","PValue","FDR")])
  dea <- aggregate(dea[,-1], by=list(symbol=dea$symbol),FUN=mean)
  rownames(dea) <- dea$symbol
  dea <- dea[,-1]
  
  lapply(dea, FUN = function(x){
    if(any(is.infinite(x))){
      w <- which(is.infinite(x))
      x[w] <- max(abs(x[setdiff(which(dea$FDR<0.5),w)]))*sign(x[w])
    }
  })
  return(dea)
}

#################################################
#' TSperm
#'
#' @param TS TargetScan DataFrame of miRNA tx targets
#' @param TP string of miRNA seed (true-positive), eg. "ACAUCAC"
#' @param genes string vector of gene names
#' @param props vector of proportions, e.g. c(.2,.3,.4,.5)
#'
#' @return TS DataFrame with permuted miRNA targets
#' 
TSperm <- function(TS, TP, genes, props){
  # get TargetScan data for each treatment miRNA family
  TS.part <- TS[TS$family==TP,]
  lapply(props, FUN=function(p){
    sn <- floor(p*nrow(TS.part))
    genes <- genes[!genes %in% TS.part$feature]
    i <- sample(seq_len(nrow(TS.part)), sn)
    j <- sample(seq_len(length(genes)), sn)
    TS.part$feature <- as.character(TS.part$feature)
    TS.part$feature[i] <- as.character(genes[j])
    TS <- rbind(TS.part, TS[TS$family!=TP,])
    return(TS)
  })
}

#################################################
#' testCombine.fisher
#'
#' @param e1 dataframe: enrichMiR test result, e.g. e.list$`miR-122`$original@res$modScore
#' @param e2 same as e1
#' @param e3 same as e1
#' @param p.colnames character vector containing column names of pvalues
#'
#' @return dataframe containing fisher-aggregated pvalues & FDR values of the supplied enrichMiR test results
#' 
testCombine.fisher <- function(e1, e2, e3, p.colnames){
  # check input
  if(!any(colnames(e1) %in% p.colnames[1]) | !any(colnames(e2) %in% p.colnames[2]) | !any(colnames(e3) %in% p.colnames[3])){
    stop("p.colnames contains invalid column names")
  }
  # miRNA families as rownames
  e <- lapply(list(e1,e2,e3), FUN=function(x){
    if(any(colnames(x) %in% "family")){
      rownames(x) <- x[,"family"]
    }
    x })
  # get common results
  fam <- intersect(intersect(rownames(e[[1]]), rownames(e[[2]])), rownames(e[[3]]))
  # combine
  p.combined <- sapply(fam, FUN=function(x) aggregation::fisher( c(e[[1]][x,p.colnames[1]], e[[2]][x,p.colnames[2]], e[[3]][x,p.colnames[3]]) )
  )
  fdr.combined <- p.adjust(p.combined, method="fdr")
  # output
  return(
    data.frame(pvalue=p.combined[names(fdr.combined[order(fdr.combined)])], 
               FDR=fdr.combined[order(fdr.combined)])
  )
}

#################################################
#' testCombine.geom
#'
#' @param e1 dataframe: enrichMiR test result, e.g. e.list$`miR-122`$original@res$modScore
#' @param e2 same as e1
#' @param e3 same as e1
#' @param p.colnames character vector containing column names of pvalues
#'
#' @return dataframe containing fisher-aggregated pvalues & FDR values of the supplied enrichMiR test results
#' 
testCombine.geom <- function(e1, e2, e3, p.colnames){
  # check input
  if(!any(colnames(e1) %in% p.colnames[1]) | !any(colnames(e2) %in% p.colnames[2]) | !any(colnames(e3) %in% p.colnames[3])){
    stop("p.colnames contains invalid column names")
  }
  # miRNA families as rownames
  e <- lapply(list(e1,e2,e3), FUN=function(x){
    if(any(colnames(x) %in% "family")){
      rownames(x) <- x[,"family"]
    }
    x })
  # get common results
  fam <- intersect(intersect(rownames(e[[1]]), rownames(e[[2]])), rownames(e[[3]]))
  # combine
  p.combined <- sapply(fam, FUN=function(x) exp( mean( c(log(e[[1]][x,p.colnames[1]]), log(e[[2]][x,p.colnames[2]]), log(e[[3]][x,p.colnames[3]])) ))
  )
  fdr.combined <- p.adjust(p.combined, method="fdr")
  # output
  return(
    data.frame(pvalue=p.combined[names(fdr.combined[order(fdr.combined)])], 
               FDR=fdr.combined[order(fdr.combined)])
  )
}

#################################################
#' doBenchmark
#'
#' @param res enrichMiR test results (e object)
#' @param TP character vector: contains the seeds (families) for a miRNA treatment
#'
#' @return a dataframe containing scores for each enrichMiR test
#'
doBenchmark <- function(res, TP){
  res <- lapply(res, FUN=function(x){
    x <- x[order(x$FDR),]
    if("family" %in% colnames(x)) row.names(x) <- x$family
    x$truth <- row.names(x) %in% TP
    x$FDR[is.na(x$FDR)] <- 1
    x
  })
  data.frame( method=names(res),
              detPPV = sapply(res, FUN=function(x) 1/which(x$truth)[1] ),
              FP.atFDR05 = sapply(res, FUN=function(x) sum(!x$truth & x$FDR<0.05)),
              log10QDiff = sapply(res, FUN=function(x){
                tp1 <- -log10(x$FDR[which(x$truth)[1]])
                fp1 <- -log10(x$FDR[which(!x$truth)[1]])
                tp1-fp1
              }),
              log10QrelIncrease = sapply(res, FUN=function(x){
                tp1 <- -log10(x$FDR[which(x$truth)[1]])
                fp1 <- -log10(x$FDR[which(!x$truth)[1]])
                (tp1-fp1)/(min(tp1,fp1))
              }),
              TP.atFDR05 = sapply(res, FUN=function(x) sum(x$truth[x$FDR<0.05]))
  )
  
}

#################################################
#' getDEAList
#'
#' @param se 
#' @param dea.df.names 
#' @param dea.names 
#' @param cleanRowNames
#'
#' @return
#' 
getDEAList <- function(se, dea.df.names, dea.names, cores=1, cleanRowNames=FALSE){
  if(length(dea.names)>1){
    dea.list <- lapply(rowData(se)[,dea.df.names], function(x) x)
  } else {
    dea.list <- list(rowData(se)[,dea.df.names])
  }
  dea.list <- lapply(dea.list, getDEA, cleanRowNames=cleanRowNames)
#  dea.list <- bplapply(dea.list, getDEA, cleanRowNames,
#                       BPPARAM = MulticoreParam(cores, progressbar = FALSE) )
  names(dea.list) <- dea.names
  return(dea.list)
}

#################################################
#' getTPbartel
#'
#' @param celltype 
#' @param dea.names 
#'
#' @return
#' @export
#'
#' @examples
getTPbartel <- function(celltype, dea.names){
  # generate TP object
  ## load detailed Bartel treatment info containing exact sequences
  pert <- read.csv("data/bartel_treatments.txt", sep = "\t")
  ## we extract HeLa & HEK data (we disregard passenger strands as long as we're working with TargetScan)
  colnames(pert) <- c("treatment","seq")
  if(celltype=="hela"){
    exclude <- pert[19:35,]
    pert <- pert[1:17,]
  } else {
    exclude <- pert[50:61,]
    pert <- pert[37:48,]
  }
  ## extract seed sequences
  pert$seed <- sapply(pert$seq, function(x) substring(x, 2, 8))
  ## get seed family for each treatment miRNA (true positive)
  mir.names <- sapply(dea.names, function(x) paste(unlist(strsplit(x,"\\.")), collapse = "-"))
  TP <- sapply(mir.names, FUN=function(x) unique(as.character(
    pert$seed[grepl(paste0(x,"\\("), pert$treatment) | grepl(paste0(x,"$"), pert$treatment)] 
  )) )
  return(TP)
}

#################################################
#' getTP
#'
#' @param tp.mirna 
#' @param TS 
#' @param dea.names 
#'
#' @return
#' @export
#'
#' @examples
getTP <- function(tp.mirna, TS, dea.names){
  # find miRNA seed sequences in TS
  TP <- sapply(tp.mirna, function(x) metadata(TS)$families[grepl(x, names(metadata(TS)$families))] )
  TP <- TP[!duplicated(TP)]
  names(TP) <- dea.names
  return(TP)
}

#################################################
#' reduceTS
#'
#' @param tp.mirna 
#' @param TS 
#'
#' @return
#' @export
#'
#' @examples
reduceTS <- function(tp.mirna, TS){
  # exclude passenger strands from consideration
  passenger <- sapply(strsplit(tp.mirna, "-"), function(x){
    guide <- rev(strsplit(x,"-"))[1]
    pass <- if(guide=="3p") "5p" else if(guide=="5p") "3p"
    paste(c(x[-3],pass),collapse = "-")
    })
  seed <- sapply(passenger, function(x)
    metadata(TS)$families[grepl(x, names(metadata(TS)$families))])
  TS <- TS[!(TS$family %in% seed),]
}



#################################################
#' getEnrichMiR
#'
#' @param dea.names 
#' @param TS 
#' @param TP 
#' @param props 
#' @param rownames.se 
#' @param tests 
#' @param dea.list 
#' @param mirexpr 
#'
#' @return
#' @export
#'
#' @examples
getEnrichMiR <- function(dea.list, TS, TP, props, rownames.se, tests, mirexpr, cores=1){
  ### do TS permutations
  set.seed(1234)
  
  if(length(grep("^ENSG|^ENSRNOG", row.names(dea.list[[1]])))>100){
    dea.list <- lapply(dea.list, FUN=function(x){
      x <- x[order(x$FDR, x$PValue),]
      g <- sapply(strsplit(row.names(x), "\\."), FUN=function(x) x[2])
      x <- x[!duplicated(g),]
      row.names(x) <- g[!duplicated(g)]
      x
    })
  }
  
  myf <- function(dea, TS){
    e <- enrichMiR(DEA=dea, TS=TS, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)
    e@TS <- data.frame() # reduce memory footprint
    e@DEA <- data.frame()
    e@res <- lapply(e@res, FUN=function(x){
      x$features <- NULL
      x
    })
    e
  }

  ### enrichMiR
  if(cores==1){
    e.list <- lapply(dea.list, FUN=function(dea){
      dea <- as.data.frame(dea)
      eres <- lapply(TSperm(TS, TP, rownames.se, props), FUN=function(ts) myf(dea, ts))
      names(eres) <- names(props)
      c( list(original=enrichMiR(DEA=dea, TS=TS, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)),
       eres )
    })
  }else{
    if(length(dea.list)==1){
      e.list <- lapply(dea.list, FUN=function(dea){
        dea <- as.data.frame(dea)
        eres <- bplapply(TSperm(TS, TP, rownames.se, props), BPPARAM=MulticoreParam(cores, progressbar=TRUE),
                         FUN=function(ts) myf(dea, ts))
        names(eres) <- names(props)
        c( list(original=enrichMiR(DEA=dea, TS=TS, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)),
           eres )
      })
    }else{
      e.list <- bplapply(dea.list, BPPARAM=MulticoreParam(cores, progressbar=TRUE), 
                         FUN=function(dea){
         dea <- as.data.frame(dea)
         eres <- lapply(TSperm(TS, TP, rownames.se, props), FUN=function(ts) myf(dea, ts))
         names(eres) <- names(props)
         c( list(original=enrichMiR(DEA=dea, TS=TS, miRNA.expression=mirexpr, cleanNames=TRUE, tests=tests)),
            eres )
      })
    }
  }

  ### combine michael, wEN and modScore results
  #### fisher aggregation
  for(i in names(e.list)){
    for(j in names(e.list[[i]])){
      comb <- testCombine.fisher(e.list[[i]][[j]]@res$michael.down, e.list[[i]][[j]]@res$wEN.down, e.list[[i]][[j]]@res$modScore, c("over.pvalue","over.pvalue","pvalue"))
      e.list[[i]][[j]]@res$combFish.1 <- comb
    }
  }
  #### geometric mean
  for(i in names(e.list)){
    for(j in names(e.list[[i]])){
      comb <- testCombine.geom(e.list[[i]][[j]]@res$michael.down, e.list[[i]][[j]]@res$wEN.down, e.list[[i]][[j]]@res$modScore, c("over.pvalue","over.pvalue","pvalue"))
      e.list[[i]][[j]]@res$combGeom.1 <- comb
    }
  }
  ### combine michael, aREAmir and KS2 results
  #### fisher aggregation
  for(i in names(e.list)){
    for(j in names(e.list[[i]])){
      comb <- testCombine.fisher(e.list[[i]][[j]]@res$michael.down, e.list[[i]][[j]]@res$aREAmir, e.list[[i]][[j]]@res$KS2, c("over.pvalue","pvalue","ks.pvalue.down"))
      e.list[[i]][[j]]@res$combFish.2 <- comb
    }
  }
  #### geometric mean
  for(i in names(e.list)){
    for(j in names(e.list[[i]])){
      comb <- testCombine.geom(e.list[[i]][[j]]@res$michael.down, e.list[[i]][[j]]@res$aREAmir, e.list[[i]][[j]]@res$KS2, c("over.pvalue","pvalue","ks.pvalue.down"))
      e.list[[i]][[j]]@res$combGeom.2 <- comb
    }
  }

  ### combine michael, aREAmir and regmirb
  #### fisher aggregation
  for(i in names(e.list)){
    for(j in names(e.list[[i]])){
      comb <- testCombine.fisher(e.list[[i]][[j]]@res$michael.down, e.list[[i]][[j]]@res$aREAmir, e.list[[i]][[j]]@res$regmirb.down, c("over.pvalue","pvalue","pvalue"))
      e.list[[i]][[j]]@res$combFish.3 <- comb
    }
  }
  #### geometric mean
  for(i in names(e.list)){
    for(j in names(e.list[[i]])){
      comb <- testCombine.geom(e.list[[i]][[j]]@res$michael.down, e.list[[i]][[j]]@res$aREAmir, e.list[[i]][[j]]@res$regmirb.down, c("over.pvalue","pvalue","pvalue"))
      e.list[[i]][[j]]@res$combGeom.3 <- comb
    }
  }

  return(e.list)
}

#################################################
#' getDF
#'
#' @param dea.names 
#' @param e.list 
#' @param TP 
#'
#' @return
#' @export
#'
#' @examples
getDF <- function(dea.names, e.list, TP){ 
  ### summarize results
  #### generate the benchmarking scores
  BM.list <- lapply(dea.names, FUN=function(x){
    lapply(names(e.list[[x]]), FUN=function(y){
      doBenchmark(e.list[[x]][[y]]@res, TP[x]) 
    })
  })
  #### naming
  names(BM.list) <- dea.names
  for(x in dea.names){
    names(BM.list[[x]]) <- names(e.list[[x]])
  }
  #### generate a results df for plotting
  BM.list2 <- lapply(BM.list, FUN=function(x) dplyr::bind_rows(x, .id = "prop.rep"))
  BM.df <- dplyr::bind_rows(BM.list2, .id="treatment")
  BM.df$prop <- unlist(lapply(strsplit(BM.df$prop.rep, "[.]"), FUN=function(x) x[1]))
  BM.df$prop <- factor(BM.df$prop, levels=unique(BM.df$prop))
  
  return(BM.df)
}
  


##########################################################################################
# main


# ## datasets to loop over
# datasets <- c("bartel","bartel.exon","jeong","jeong.exon","ratPolyA","ratPolyA.exon")
## cores
cores <- 6
## set number of replicates per permutation proportion
nrep <- 3
names(i) <- i <- 1:nrep
props <- unlist(lapply( c("20"=0.2, "35"=0.35, "50"=0.5), 
                        FUN=function(x) lapply(i, FUN=function(y) x)))
## which enrichMiR tests to run (NULL for all)
#tests <- c("overlap","michael","modscore","modsites","areamir","regmirb")



TPs <- c( 
  let.7a = "GAGGUAG", lsy.6 = "UUUGUAU", miR.1 = "GGAAUGU", miR.124 = "AAGGCAC", 
  miR.137 = "UAUUGCU", miR.139 = "CUACAGU", miR.143 = "GAGAUGA", 
  miR.144 = "ACAGUAU", miR.153 = "UGCAUAG", miR.155 = "UAAUGCU", 
  miR.182 = "UUGGCAA", miR.199a = "CCAGUGU", miR.204 = "UCCCUUU", 
  miR.205 = "CCUUCAU", miR.216b = "AAUCUCU", miR.223 = "GUCAGUU", 
  miR.7 = "GGAAGAC", miR.122 = "GGAGUGU", miR.133 = "UGGUCCC", 
  miR.138 = "GCUGGUG", miR.145 = "UCCAGUU", miR.184 = "GGACGGA", 
  miR.190a = "GAUAUGU", miR.200b = "AAUACUG", miR.216a = "AAUCUCA", 
  miR.217 = "ACUGCAU", miR.219a = "GAUUGUC", miR.375 = "UUGUUCG", 
  miR.451a = "AACCGUU", "DKOvWT"="GCUACAU", "DKO"="GCUACAU", 
  "miR.138vNeg"="GCUGGUG", "miR.499vNeg"="UAAGACU", 
  "miR.138"="GCUGGUG", "miR.499"="UAAGACU")

allds=c("bartel.hela.DEA.SE.rds",
       "bartel.hek.DEA.SE.rds",
       "bartel.hela.exon.DEA.SE.rds",
       "bartel.hek.exon.DEA.SE.rds",
       "jeong.DEA.SE.rds",
       "jeong.exon.DEA.SE.rds",
       "ratPolyA.DEA.SE.rds",
       "ratPolyA.exon.DEA.SE.rds")


# miRNAs
getMirexpr <- function(tissue){
  data("microRNAome")
  counts <- rowMedians(cpm(calcNormFactors(DGEList(assay(microRNAome[,microRNAome$cell_tissue==tissue])))))
  names(counts) <- rownames(microRNAome)
  return(counts[counts>=10])
}

hela <- getMirexpr("hela")
hek <- getMirexpr("hek293")
ratmirna <- read.csv("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/DIV20_Neuron_HC_miRNA.csv")
rat <- 2^(ratmirna$logCPM) 
names(rat) <- ratmirna$Gene.name
rat <- rat[rat>=10]

mirexpr <- list(hela=hela, 
                hek=hek, 
                rat=rat)


allParams <- c(
  lapply(allds[1:6], FUN=function(x) 
  list(ds=x, ts="TargetScan_human.rds", miRNAs=NULL) ),
  lapply(allds[7:8], FUN=function(x) 
    list(ds=x, ts="TargetScan_rat.rds", miRNAs=NULL) )
)

mirParams <- c(
  lapply(allds[c(1,3)], FUN=function(x) 
    list(ds=x, ts="TargetScan_human.rds", miRNAs=mirexpr$hela) ),
  lapply(allds[c(2,4)], FUN=function(x) 
    list(ds=x, ts="TargetScan_human.rds", miRNAs=mirexpr$hek) ),
  lapply(allds[7:8], FUN=function(x) 
    list(ds=x, ts="TargetScan_rat.rds", miRNAs=mirexpr$rat) )
)


doRun <- function(p, cores=6){
  se <- readRDS(file.path("data", p$ds))
  dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))])
  dea.names <- gsub("DEA\\.","", dea.df.names)
  dea.names <- gsub("spliced\\.","", dea.names)
  dea.names <- gsub("^\\.","", dea.names)
  w <- which(dea.names %in% names(TPs))
  if(length(w)==0) stop("Unknown DEAs (", paste(dea.names, collapse=", "), ")")
  dea.names <- dea.names[w]
  dea.df.names <- dea.df.names[w]
  ### get list of all DEAs
  dea.list <- getDEAList(se, dea.df.names, dea.names)
  ### specify expressed miRNA (celltype-specific)
  mirexpr <- p$miRNAs
  ### load TS object
  TS <- readRDS( file.path("data",p$ts) )
  TS$sites <- round(TS$sites)
  ### get true positive miRNA treatments
  TP <- TPs[dea.names]
  ### enrichMiR
  eres <- getEnrichMiR( dea.list=dea.list, TS=TS, TP=TP, props = props, tests=tests,
                        rownames.se=rownames(se), mirexpr=p$miRNAs, cores=cores )
  ### summarize enrichMiR results in a dataframe
  bench <- getDF(dea.names, eres, TP)
  ### output
  if(is.null(p$miRNAs)){
    saveRDS(eres, file=file.path("results", gsub("DEA\\.SE","enrichMiR",p$ds)))
    saveRDS(bench, file=file.path("results", gsub("DEA\\.SE","benchmark",p$ds)))
  } else {
    saveRDS(eres, file=file.path("results", gsub("DEA\\.SE","mirexpr.enrichMiR",p$ds)))
    saveRDS(bench, file=file.path("results", gsub("DEA\\.SE","mirexpr.benchmark",p$ds)))
  }
  rm(eres)
  gc(verbose = FALSE, full = TRUE)
  
}


