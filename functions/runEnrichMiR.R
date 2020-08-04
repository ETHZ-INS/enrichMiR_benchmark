#--------------------------------------------------------------------------------
#' getDEAList
#'
#' @param p list at least containing elements 'dea' and 'ann'
#'
#' @return p list with loaded DEAs & annotation and miRNA library in $families
#' 
loadAll <- function(p){
  
  # DEA
  
  if(is.character(p$dea)){
    library(SummarizedExperiment)
    ## load file
    dea <- readRDS(p$dea)
    ## if the file is an SE extract the DEAs
    if(!is.list(dea)){
      se <- dea
      dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))])
      dea.names <- gsub("DEA\\.","", dea.df.names)
      dea.names <- gsub("spliced\\.","", dea.names)
      dea.names <- gsub("^\\.","", dea.names)
      w <- which(dea.names %in% names(p$tp))
      if(length(w)==0) stop("Unknown DEAs (", paste(dea.names, collapse=", "), ")")
      dea.names <- dea.names[w]
      dea.df.names <- dea.df.names[w]
      
      if(length(dea.names)>1){
        dea.list <- lapply(rowData(se)[,dea.df.names], function(x) x)
      } else {
        dea.list <- list(rowData(se)[,dea.df.names])
      }
      names(dea.list) <- dea.names
      p$dea <- dea.list
      ## if the file is a list of DEAs just append it to p
    } else {
      p$dea <- dea
    }
  }
  
  # Annotation
  
  if(is.character(p$ann)){
    ## load Annotation
    ann <- readRDS(p$ann)
    ann$sites <- round(ann$sites)
    p$ann <- ann
  }
  if(!is.null(metadata(p$ann)$families)){
    ## add miRNA library to p
    p$families <- metadata(p$ann)$families
    metadata(p$ann)$families <- NULL
  }
  
  return(p)
}


#--------------------------------------------------------------------------------
#' cleanDEA
#'
#' @param p list at least containing element 'dea'
#'
#' @return p list with DEA rownames that should correspond to annotation features
#'
cleanDEA <- function(p){
  
  p$dea <- lapply(p$dea, function(dea){
    
    # find correct column names & order
    fdr.n <- if("padj" %in% colnames(dea)) "padj" else if("FDR" %in% colnames(dea)) "FDR"
    pval.n <- if("pvalue" %in% colnames(dea)) "pvalue" else if("PValue" %in% colnames(dea)) "PValue"
    dea <- dea[order(dea[,fdr.n], dea[,pval.n]),]
    # when using gene-level annotations: clean up rownames & aggregate tx
    if(length(grep("^ENSG|^ENSRNOG", row.names(dea))) > 100){
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
    
    dea
  })
  
  return(p)
}


#--------------------------------------------------------------------------------
#' TSPerm
#'
#' @param p list at least containing elements 'dea', 'ann', 'props'
#'
#' @return p list containing permuted annotations based on 'props'
#' 
TSPerm <- function(p){
  
  TS <- p$ann
  genes <- rownames(p$dea[[1]])
  props <- p$props
  
  # get TargetScan data for each treatment miRNA family
  if(class(TS) %in% c("DFrame","data.frame")){
    TS.list <- lapply(p$tp, function(TP){
      TS.part <- TS[TS$family == TP,]
      lapply(props, FUN=function(p){
        sn <- floor(p*nrow(TS.part))
        genes <- genes[!genes %in% TS.part$feature]
        i <- sample(seq_len(nrow(TS.part)), sn)
        j <- sample(seq_len(length(genes)), sn)
        TS.part$feature <- as.character(TS.part$feature)
        TS.part$feature[i] <- as.character(genes[j])
        TS <- rbind(TS.part, TS[TS$family!=TP,])
        TS
      })
    })
    
  } else if (class(TS) %in% "GRanges"){
    
    #seeds <- sapply(levels(tar.ann$seed),function(x) metadata(TS)$families[x])
    #names(seeds) <- levels(tar.ann$seed)
  }
  
  p$ann <- c(original=TS, TS.list)
  return(p)
}


#--------------------------------------------------------------------------------
#' getEnrichMiR
#'
#' @param p list at least containing elements 'dea', 'ann', 'tp', 'tests', 'mirnas', 'perm', 'props' and 'families'
#' @param cores integer; cores for bplapply()
#'
#' @return list of enrichMiR results
#'
getEnrichMiR <- function(p, cores=1){
  
  # enrichMiR function
  
  callEnrichMiR <- function(dea, TS){
    e <- enrichMiR(DEA=dea, TS=TS, miRNA.expression=p$mirnas, cleanNames=FALSE, tests=p$tests, families=p$families)
    # reduce memory footprint
    e@TS <- data.frame() 
    e@DEA <- data.frame()
    e@res <- lapply(e@res, FUN=function(x){
      x$features <- NULL
      x
    })
    return(e)
  }
  
  # main 
  
  ## define miRNA library (either organism-specific or tissue-specific if mirexpr supplied)
  if(!is.null(p$mirnas)){
    p$families <- p$families[intersect(names(p$families), names(p$mirnas))]
  }
  
  ## enrichMiR with annotation permutations
  if(p$perm){
    if(length(p$dea)==1){
      e.list <- lapply(names(p$tp), function(tp){
        # enrichMiR for all annotation props
        eres <- bplapply(p$ann[[tp]], BPPARAM=MulticoreParam(cores, progressbar=TRUE), 
                         function(ann) callEnrichMiR(as.data.frame(p$dea[[tp]]), ann) )
        names(eres) <- names(p$props)
        # combine with enrichMiR for original annotation 
        eres <- c(list(
          original=enrichMiR(DEA=as.data.frame(p$dea[[tp]]), TS=p$ann$original, miRNA.expression=p$mirnas, 
                             cleanNames=FALSE, tests=p$tests, families=p$families) ), eres)
      })
    } else {
      e.list <- bplapply(names(p$tp), BPPARAM=MulticoreParam(cores, progressbar=TRUE), function(tp){
        # enrichMiR for all annotation props
        eres <- lapply(p$ann[[tp]], function(ann) callEnrichMiR(as.data.frame(p$dea[[tp]]), ann) )
        names(eres) <- names(p$props)
        # combine with enrichMiR for original annotation 
        eres <- c(list(
          original=enrichMiR(DEA=as.data.frame(p$dea[[tp]]), TS=p$ann$original, miRNA.expression=p$mirnas, 
                             cleanNames=FALSE, tests=p$tests, families=p$families) ), eres)
      })
    }
    names(e.list) <- names(p$tp)
    
    ## enrichMiR without annotation permutations
  } else {
    if(length(p$dea)==1){
      e.list <- lapply(names(p$tp), function(tp){
        enrichMiR(DEA=as.data.frame(p$dea[[tp]]), TS=p$ann, miRNA.expression=p$mirnas, 
                  cleanNames=FALSE, tests=p$tests, families=p$families) 
      })
    } else {
      e.list <- bplapply(names(p$tp), BPPARAM=MulticoreParam(cores, progressbar=TRUE), function(tp){
        enrichMiR(DEA=as.data.frame(p$dea[[tp]]), TS=p$ann, miRNA.expression=p$mirnas, 
                  cleanNames=FALSE, tests=p$tests, families=p$families)
      })
    }
    names(e.list) <- names(p$tp)
  }
  
  return(e.list)
}


#--------------------------------------------------------------------------------
#' testCombine.fisher
#'
#' @param e DataFrame; contains enrichMiR test results for a single treatment
#' @param pval.n named string vector; contains the column name harbouring pvalues for each test that should be aggregated
#'
#' @return a single DataFrame containing the columns pvalues & FDR of the fisher aggregation
#' 
testCombine.fisher <- function(e, pval.n){
  
  # get common families
  fam <- lapply(e, function(x) rownames(x))
  fam <- Reduce(intersect, fam)
  
  # get a matrix of p-values for each of the tests we want to combine
  p <- sapply(1:length(e), function(i) sapply(fam, function(j) e[[i]][j,pval.n[[i]]] ))
  # combine the p-values for each family; creates vector
  p.comb <- sapply(fam, function(i) aggregation::fisher(p[i,]) )
  # calculate FDR values
  fdr.comb <- p.adjust(p.comb, method="fdr")
  
  return( data.frame(pvalue=p.comb[names(fdr.comb[order(fdr.comb)])], 
                     FDR=fdr.comb[order(fdr.comb)]) 
  )
}


#--------------------------------------------------------------------------------
#' testCombine.geom
#'
#' @param e DataFrame; contains enrichMiR test results for a single treatment
#' @param pval.n named string vector; contains the column name harbouring pvalues for each test that should be aggregated
#'
#' @return a single DataFrame containing the columns pvalues & FDR of the geometric mean aggregation
#' 
testCombine.geom <- function(e, pval.n){
  
  # get common families
  fam <- lapply(e, function(x) rownames(x))
  fam <- Reduce(intersect, fam)
  
  # get a matrix of p-values for each of the tests we want to combine
  p <- sapply(1:length(e), function(i) sapply(fam, function(j) e[[i]][j,pval.n[[i]]] ))
  # combine the p-values for each family; creates vector
  p.comb <- sapply(fam, function(i) exp(mean(log( p[i,] ))))
  # calculate FDR values
  fdr.comb <- p.adjust(p.comb, method="fdr")
  
  return( data.frame(pvalue=p.comb[names(fdr.comb[order(fdr.comb)])], 
                     FDR=fdr.comb[order(fdr.comb)])
  )
}


#--------------------------------------------------------------------------------
#' combineTests
#'
#' @param e.list list of DataFrames; contains results of enrichMiR tests
#' @param tests list of named string vectors; each string is the name of a DataFrame in e.list that should be used for aggregation; if no list is supplied 3 default test combinations will be aggregated
#' @param fisher.agg binary; TRUE for fisher aggregation
#' @param geom.agg binary; TRUE for geometric mean aggregation
#'
#' @return the initial e.list with added aggregated enrichMiR test results; each aggregation consists of pvalue & FDR columns
#' 
combineTests <- function(e.list, tests=NULL, fisher.agg=TRUE, geom.agg=TRUE){
  
  # default enrichMiR tests to combine 
  if(is.null(tests)){
    tests <- list(
      comb1=c("siteMir.down","wEN.down","modScore"),
      comb2=c("siteMir.down","aREAmir","KS2"),
      comb3=c("siteMir.down","aREAmir","regmirb.down")
    )
  }
  
  # possible names of columns containing p-Values of interest
  pval.cols <- c("over.pvalue","pvalue","ks.pvalue.down")
  
  # actual names of p-Value columns based on test outputs
  pval.n <- lapply(tests, function(set){
    unlist( sapply(set, function(i){
      unlist( lapply(pval.cols, function(x){
        e <- e.list[[1]]@res[[i]]
        colnames(e)[colnames(e)==x]
      }) )
    }) )
  })
  
  # run fisher aggregation
  if(fisher.agg) comb.f <- lapply(e.list, function(e) lapply(names(tests), function(t) DataFrame( testCombine.fisher(e@res[tests[[t]]], pval.n[[t]]) )))
  # run geometric mean aggregation
  if(geom.agg) comb.g <- lapply(e.list, function(e) lapply(names(tests), function(t) testCombine.geom(e@res[tests[[t]]], pval.n[[t]]) ))
  
  
  # transform comb outputs into DFrames; append to e.list via custom names (names(tests))
  
  # add names to aggregated tests
  for(i in names(e.list)){
    if(fisher.agg) names(comb.f[[i]]) <- paste0(names(tests), ".fish")
    if(geom.agg) names(comb.g[[i]]) <- paste0(names(tests), ".geom")
  }
  
  # add aggregation results to e.list
  for(i in names(e.list)){
    if(fisher.agg) e.list[[i]]@res[names(comb.f[[i]])] <- comb.f[[i]]
    if(geom.agg) e.list[[i]]@res[names(comb.g[[i]])] <- comb.g[[i]]
  }
  
  return(e.list)
}

