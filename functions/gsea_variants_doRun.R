####################################################################################################
####################################################################################################
# FUNCTIONS: common

devtools::load_all("../../tgermade_test/master_19_20/enrichMiR/enrichMiR/")

#' cleanDEA
#'
#' @param dea 
#'
#' @return
#' @export
#'
#' @examples
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

#' loadAll
#'
#' @param dea 
#' @param tp 
#'
#' @return
#' @export
#'
#' @examples
loadAll <- function(dea,tp){
  if(is.character(dea)){
    library(SummarizedExperiment)
    ## load file
    dea <- readRDS(dea)
    ## if the file is an SE extract the DEAs
    if(!is.list(dea)){
      se <- dea
      dea.df.names <- colnames(rowData(se)[,grepl("DEA",colnames(rowData(se)))])
      dea.names <- gsub("DEA\\.","", dea.df.names)
      dea.names <- gsub("spliced\\.","", dea.names)
      dea.names <- gsub("^\\.","", dea.names)
      w <- which(dea.names %in% names(tp))
      if(length(w)==0) stop("Unknown DEAs (", paste(dea.names, collapse=", "), ")")
      dea.names <- dea.names[w]
      dea.df.names <- dea.df.names[w]
      
      if(length(dea.names)>1){
        dea.list <- lapply(rowData(se)[,dea.df.names], function(x) x)
      } else {
        dea.list <- list(rowData(se)[,dea.df.names])
      }
      names(dea.list) <- dea.names
      dea <- dea.list
    } else {
      dea <- dea
    }
  }
  return(dea)
}


####################################################################################################
####################################################################################################
# FUNCTIONS: To run for combined tests (results from other tests + 2nd regmir test: different variations)

#' gsea.new
#'
#' @param signal A named logical vector indicating which features are differentially-expressed
#' @param sets A data.frame with at least the following columns: 'set', 'feature'.
#' @param maxSize The maximum number of elements in a set to be tested (default 500). If the test takes too long to run, consider setting this.
#' @param nperm The number of permutations, default 2000. The more permutations, the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @importFrom fgsea fgsea
#' @export
gsea.new <- function(e, dea, sets, maxSize=300, nperm=2000, KO=TRUE, useFC=TRUE, sig.binary=TRUE, pure=FALSE, ...){
  # generate signal
  if(sig.binary){
    binary.signatures <- list(
      down=.dea2binary(dea, th=0.05, th.alfc=0, restrictSign=-1),
      up=.dea2binary(dea, th=0.05, th.alfc=0, restrictSign=1)
    )
    if(KO){
      signal <- binary.signatures$up
    } else {
      signal <- binary.signatures$down
    }
  } else {
    if(useFC){
      signal <- dea$logFC
    } else {
      dea$FDR[dea$FDR==0] <- min(dea$FDR[dea$FDR>0])
      signal <- -log10(dea$FDR)*sign(dea$logFC)
    }
    names(signal) <- rownames(dea)
  }
  # clean e & restrict sets
  if(!pure){
    e <- as.data.frame(e)
    e[is.na(e)] <- 1
    # generate co OLD
    # if( sum(e$FDR<=.05,na.rm=TRUE) > 30 ){
    #   co <- rownames(e[e$FDR<=.05,])
    # } else if( sum(e$FDR<=.2) > 30 ){
    #   co <- rownames(e[e$FDR<=.2,])
    # } else if( nrow(e)>=30 ) {
    #   co <- rownames(e[1:30,])
    # } else {
    #   co <- rownames(e[1:nrow(e),])
    # }
    # generate co NEW
    if( nrow(e) >= 30 ){
      co <- head(rownames(e),n=max(30,ceiling(nrow(e)/10)))
    } else {
      co <- rownames(e)
    }
    sets <- as.data.frame(sets[sets$set %in% co,])
    if(is.factor(sets$set)) sets$set <- droplevels(sets$set)
  }
  
  sets <- lapply(split(sets$feature,sets$set), tested=names(signal), 
                 FUN=function(x,tested){ intersect(unique(x),tested) })
  res <- fgsea::fgsea(sets, signal, nperm, minSize=4, maxSize=maxSize, ...)
  #print("gsea run successful.")
  if(nrow(res)==0) return("No GSEA output.")
  res <- as.data.frame(res)
  res <- res[order(res$padj,res$pval),]
  colnames(res)[1:5] <- c("family","pvalue","FDR","ES","normalizedEnrichment")
  colnames(res)[8] <- "features"
  row.names(res) <- res[,1]
  return(res[,-1])
}


#' gseaRun
#'
#' @param p 
#'
#' @return
#' @export
#'
#' @examples
gseaRun <- function(p, pure=FALSE){
  dea <- p$dea
  dea <- loadAll(dea, TPs)
  dea <- lapply(dea, cleanDEA)
  e.list <- readRDS(p$e)
  TS <- readRDS(p$ts)
  if(all(c("family", "feature") %in% colnames(TS))){
    colnames(TS) <- gsub("^family$","set", colnames(TS))
  }
  TS <- DataFrame(TS)
  ko <- p$ko
  n <- p$n
    
  print(paste("gsea pure:", pure))
  
  for(bin in c(T,F)){
    S=ifelse(bin,"bin","cont")
    if(bin) use.fc=F else use.fc=c(T,F)
    
    for(fc in use.fc){
      print(paste("signal type:", toupper(S)))
      if(!bin){
        FC=ifelse(fc,"fc","fdr")
        print(paste("signal feature:", toupper(FC)))
      }
      
      if(pure){
        res <- lapply(names(dea), function(i) 
          gsea.new(e=NULL, dea[[i]], TS, KO=ko, sig.binary=bin, useFC=fc, pure=pure, BPPARAM=MulticoreParam(8))
        )
        dir <- "../results_gsea/gsea_pure/"
      } else {
        res <- lapply(names(dea), function(i) lapply(e.list[[i]]$original@res, function(e)
          gsea.new(e, dea[[i]], TS, KO=ko, sig.binary=bin, useFC=fc, pure=pure, BPPARAM=MulticoreParam(8))
        ))
        dir <- "../results_gsea/"
      }
      names(res) <- names(dea)
      
      if(bin){
        saveRDS(res, paste0(dir, n, ".gsea.", S,".rds") )
        print( paste(paste0(n, ".gsea.", S,".rds"),"saved.") )
      } else {
        saveRDS(res, paste0(dir, n, ".gsea.", S,"_",FC,".rds") )
        print( paste(paste0(n, ".gsea.", S,"_",FC,".rds"),"saved.") )
      }
    }
  }
}


####################################################################################################
####################################################################################################
# MAIN

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
  "miR.138vNeg"="GCUGGUG", "miR.499vNeg"="UAAGACU", "miR.499"="UAAGACU", 
  "218DKOvWT"="UGUGCUU", "218DKO"="UGUGCUU",
  "138DKOvWT"="GCUGGUG", "138DKO"="GCUGGUG"
)

# load
n <- c("bartel.hek","bartel.hela","amin","ratPolyA","cherone.d0","cherone.d1","jeong")
dea.files <- paste0("../data/",n,".DEA.SE.rds")
e.files <- paste0("../results/",n,".enrichMiR.rds")
ts.files <- list.files("../data", pattern="TargetScan_", full.names=TRUE)
names(ts.files) <- c("hsa","mmu","rno")

params <- lapply(1:length(n), function(i){
  if( grepl("bartel|jeong",n[i]) ){
    org="hsa"
  } else if( grepl("ratPolyA",n[i]) ){
    org="rno"
  } else org="mmu"
  if( grepl("bartel|ratPolyA",n[i]) ){
    ko=FALSE
  } else ko=TRUE
  list(n=n[i], dea=dea.files[i], e=e.files[i], ts=ts.files[org], ko=ko)
})
names(params) <- n

# run: different regmir variations based on other tests
#lapply(params, gseaRun)

# run: different regmir variations
lapply(params, gseaRun, pure=TRUE)
