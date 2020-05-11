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
})

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
getDEAList <- function(se, dea.df.names, dea.names, cores, cleanRowNames=FALSE){
  if(length(dea.names)>1){
    dea.list <- lapply(rowData(se)[,dea.df.names], function(x) x)
  } else {
    dea.list <- list(rowData(se)[,dea.df.names])
  }
  dea.list <- bplapply(dea.list, getDEA, cleanRowNames,
                       BPPARAM = MulticoreParam(cores, progressbar = FALSE) )
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
getEnrichMiR <- function(dea.names, TS, TP, props, rownames.se, tests, dea.list, mirexpr){
  ### do TS permutations
  set.seed(1234)
  TS.list <- bplapply(dea.names, FUN=function(i){
    TSperm(TS, TP[i], rownames.se, props)
  }, BPPARAM = MulticoreParam(cores, progressbar = FALSE) )
  
  ### naming
  names(TS.list) <- dea.names
  for(i in dea.names){
    names(TS.list[[i]]) <- names(props)
  }
  ### add original TS to permutated ones
  for(i in dea.names){
    TS.list[[i]][["original"]] <- TS
  }
  ### place original at first place
  for(i in dea.names){
    TS.list[[i]] <- TS.list[[i]][c("original", names(props))]
  }
  ### add original to proportion vector
  props.all <- c(original="original", props)
  
  ### enrichMiR
  #### serial run
  e.list <- lapply(dea.names, FUN=function(i){
    lapply(names(props.all), FUN=function(j){
      enrichMiR(DEA=as.data.frame(dea.list[[i]]), TS=TS.list[[i]][[j]], miRNA.expression=mirexpr, 
                cleanNames=TRUE, tests=tests)
    })
  })
  #### naming
  names(e.list) <- dea.names
  for(i in dea.names){
    names(e.list[[i]]) <- names(props.all)
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


## datasets to loop over
datasets <- c("bartel","bartel.exon","jeong","jeong.exon","ratPolyA","ratPolyA.exon")
## cores
cores <- 8
## set number of replicates per permutation proportion
nrep <- 3
names(i) <- i <- 1:nrep
props <- unlist(lapply( c("20"=0.2, "30"=0.3, "40"=0.4, "50"=0.5), 
                        FUN=function(x) lapply(i, FUN=function(y) x)))
## which enrichMiR tests to run (NULL for all)
tests <- NULL
## initiate enrichMiR results
e.list <- list()
## initiate dataframe summary of results
df.list <- list()
### load enrichMiR package
devtools::load_all("/mnt/schratt/tgermade_test/master_19_20/enrichMiR/enrichMiR/")


## enrichMiR loop
for( i in datasets ){
  tryCatch({
    
    ### Bartel dataset
    ##################
    if(i=="bartel"){
      
      celltype <- c("hela", "hek")
      
      for(j in celltype){
        if(j=="hela"){
          se <- readRDS("data/bartel.hela.DEA.SE.rds")
          dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-c(1,3)])
          dea.names <- gsub("DEA.","", dea.df.names)
          #### get list of all DEAs
          dea.list <- getDEAList(se, dea.df.names,dea.names, cores)
          #### specify expressed miRNA (celltype-specific)
          mirexpr <- NULL
          #### load TS object
          TS <- readRDS( "data/TargetScan_human.rds")
          #### get true positive miRNA treatments
          TP <- getTPbartel(j, dea.names)

        } else if(j=="hek"){
          se <- readRDS("data/bartel.hek.DEA.SE.rds")
          dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-1])
          dea.names <- gsub("DEA.","", dea.df.names)
          #### get list of all DEAs
          dea.list <- getDEAList(se, dea.df.names,dea.names, cores)
          #### specify expressed miRNA (celltype-specific)
          mirexpr <- NULL
          #### load TS object
          TS <- readRDS( "data/TargetScan_human.rds")
          #### get true positive miRNA treatments
          TP <- getTPbartel(j, dea.names)
        }
        ### enrichMiR
        e.list[[i]][[j]] <- getEnrichMiR(dea.names, TS, TP, props, rownames(se), tests, dea.list, mirexpr)
        ### summarize enrichMiR results in a dataframe
        df.list[[i]][[j]] <- getDF(dea.names, e.list[[i]][[j]], TP)
          
      }

    ### Bartel exon dataset
    #######################
    } else if(i=="bartel.exon"){
      
      celltype <- c("hela", "hek")
      
      for(j in celltype){
        if(j=="hela"){
          se <- readRDS("data/bartel.hela.exon.DEA.SE.rds")
          dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-c(1,3)])
          dea.names <- gsub("DEA.spliced.","", dea.df.names)
          #### get list of all DEAs
          dea.list <- getDEAList(se, dea.df.names, dea.names, cores, cleanRowNames=TRUE)
          #### specify expressed miRNA (celltype-specific)
          mirexpr <- NULL
          #### load TS object
          TS <- readRDS( "data/TargetScan_human.rds")
          #### get true positive miRNA treatments
          TP <- getTPbartel(j, dea.names)
          
        } else if(j=="hek"){
          se <- readRDS("data/bartel.hek.exon.DEA.SE.rds")
          dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-1])
          dea.names <- gsub("DEA.spliced.","", dea.df.names)
          #### get list of all DEAs
          dea.list <- getDEAList(se, dea.df.names, dea.names, cores, cleanRowNames=TRUE)
          #### specify expressed miRNA (celltype-specific)
          mirexpr <- NULL
          #### load TS object
          TS <- readRDS( "data/TargetScan_human.rds")
          #### get true positive miRNA treatments
          TP <- getTPbartel(j, dea.names)
        }
        ### enrichMiR
        e.list[[i]][[j]] <- getEnrichMiR(dea.names, TS, TP, props, rownames(se), tests, dea.list, mirexpr)
        ### summarize enrichMiR results in a dataframe
        df.list[[i]][[j]] <- getDF(dea.names, e.list[[i]][[j]], TP)
      }

    ### Jeong dataset
    #################
    } else if(i=="jeong"){
      
      se <- readRDS("data/jeong.DEA.SE.rds")
      dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-1])
      dea.names <- "miR.221.222"
      #### get list of all DEAs
      dea.list <- getDEAList(se, dea.df.names,dea.names, cores, cleanRowNames=TRUE)
      #### specify expressed miRNA (celltype-specific)
      mirexpr <- NULL
      #### load TS object
      tp.mirna <- c("miR-221-3p","miR-222-3p")
      TS <- readRDS("data/TargetScan_human.rds")
      TS <- reduceTS(tp.mirna, TS)
      #### get true positive miRNA treatments
      TP <- getTP(tp.mirna, TS, dea.names)
      
      ### enrichMiR
      e.list[[i]] <- getEnrichMiR(dea.names, TS, TP, props, rownames(se), tests, dea.list, mirexpr)
      ### summarize enrichMiR results in a dataframe
      df.list[[i]] <- getDF(dea.names, e.list[[i]], TP)
   
    ### Jeong exon dataset
    ######################
    } else if(i=="jeong.exon"){
      
      se <- readRDS("data/jeong.exon.DEA.SE.rds")
      dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-1])
      dea.names <- "miR.221.222"
      #### get list of all DEAs
      dea.list <- getDEAList(se, dea.df.names,dea.names, cores, cleanRowNames=TRUE)
      #### specify expressed miRNA (celltype-specific)
      mirexpr <- NULL
      #### load TS object
      tp.mirna <- c("miR-221-3p","miR-222-3p")
      TS <- readRDS("data/TargetScan_human.rds")
      TS <- reduceTS(tp.mirna, TS)
      #### get true positive miRNA treatments
      TP <- getTP(tp.mirna, TS, dea.names)
      
      ### enrichMiR
      e.list[[i]] <- getEnrichMiR(dea.names, TS, TP, props, rownames(se), tests, dea.list, mirexpr)
      ### summarize enrichMiR results in a dataframe
      df.list[[i]] <- getDF(dea.names, e.list[[i]], TP)
      
    ### ratPolyA dataset
    ####################
    } else if(i=="ratPolyA"){
      
      se <- readRDS("data/ratPolyA.DEA.SE.rds")
      dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-1])
      dea.names <- levels(se$miRNA)[-1]
      #### get list of all DEAs
      dea.list <- getDEAList(se, dea.df.names,dea.names, cores, cleanRowNames=TRUE)
      #### specify expressed miRNA (celltype-specific)
      mirexpr <- NULL
      tp.mirna <- c("miR-138-5p","miR-499-5p")
      TS <- readRDS("data/TargetScan_rat.rds")
      TS <- reduceTS(tp.mirna, TS)
      #### get true positive miRNA treatments
      TP <- getTP(tp.mirna, TS, dea.names)
      
      ### enrichMiR
      e.list[[i]] <- getEnrichMiR(dea.names, TS, TP, props, rownames(se), tests, dea.list, mirexpr)
      ### summarize enrichMiR results in a dataframe
      df.list[[i]] <- getDF(dea.names, e.list[[i]], TP)
      
      
    ### ratPolyA exon dataset
    #########################
    } else if(i=="ratPolyA.exon"){
      
      se <- readRDS("data/ratPolyA.exon.DEA.SE.rds")
      dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))][-1])
      dea.names <- levels(se$miRNA)[-1]
      #### get list of all DEAs
      dea.list <- getDEAList(se, dea.df.names,dea.names, cores, cleanRowNames=TRUE)
      #### specify expressed miRNA (celltype-specific)
      mirexpr <- NULL
      tp.mirna <- c("miR-138-5p","miR-499-5p")
      TS <- readRDS("data/TargetScan_rat.rds")
      TS <- reduceTS(tp.mirna, TS)
      #### get true positive miRNA treatments
      TP <- getTP(tp.mirna, TS, dea.names)
      
      ### enrichMiR
      e.list[[i]] <- getEnrichMiR(dea.names, TS, TP, props, rownames(se), tests, dea.list, mirexpr)
      ### summarize enrichMiR results in a dataframe
      df.list[[i]] <- getDF(dea.names, e.list[[i]], TP)
    }
    
    ## save output
    saveRDS(e.list,"enrichMiR.list.rds")
    saveRDS(df.list,"benchmark_df.list.rds")
    
  }, error=function(e){
    message("Error in ", i)
    print(e)
  })
}

