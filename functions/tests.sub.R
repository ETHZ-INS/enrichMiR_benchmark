cleanDEA <- function(DEA){
  DEA <- DEA[which(!is.na(DEA$logFC)),]
  w <- which(is.infinite(DEA$logFC))
  if(length(w)>0) DEA$logFC[w] <- max(abs(DEA$logFC[-w]))*sign(DEA$logFC[w])
  if(any(DEA$FDR==0)){
    if(any(DEA$FDR<.Machine$double.xmin)){
      DEA$FDR[DEA$FDR==0] <- min(DEA$FDR[DEA$FDR!=0])
    } else {
      DEA$FDR[DEA$FDR==0] <- .Machine$double.xmin
    }
  }
  return(DEA)
}


.overlap.prob <- function (set1, set2, universe, lower = F){
  set1 <- as.character(set1)
  set2 <- as.character(set2)
  if (class(universe) == "character") {
    set1 <- intersect(set1, universe)
    set2 <- intersect(set2, universe)
    universe <- length(unique(universe))
  }
  set1 <- unique(set1)
  set2 <- unique(set2)
  ov <- sum(set1 %in% set2)
  phyper(max(0, ov - 1), length(set1), universe - length(set1), 
         length(set2), lower.tail = lower)
}



.censorScore <- function(x){
  x <- 0.1-x
  w <- which(x>0.9)
  x[w] <- 0.9+ecdf(x[w])(x[w])/10
  x
}

TS2regulon <- function(x, likelihood="score"){
  x <- as.data.frame(x)
  if(likelihood=="score"){
    x$likelihood <- .censorScore(x[[likelihood]])
  }else{
    x$likelihood <- x[[likelihood]]
  }
  lapply(split(x,x$family, drop=TRUE), FUN=function(x){
    y <- list(  tfmode=rep(-1,nrow(x)),
                likelihood=x$likelihood )
    lapply(y, FUN=function(a){
      names(a) <- x$feature
      a
    })	
  })
}


#' aREAmir
#'
#' analytic Rank-based Enrichment Analysis using a conversion of targetScan 
#' scores as weights.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of elements in a set to be tested (default 5).
#' @param maxSize The maximum number of elements in a set to be tested (default 500). If the test takes too long to run, consider setting this.
#' @param fdr.thres The FDR threshold below which genes are considered; default 0.2.
#' @param nperm The number of permutations, default 2000. The more permutations, the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @export
aREAmir <- function(dea, TS, minSize=5, pleiotropy=FALSE, useFC=TRUE){
  print("aREAmir...")
  dea <- cleanDEA(dea)
  if(useFC){
    sig <- dea$logFC
  } else {
    sig <- -log10(dea$FDR)*sign(dea$logFC)
  }
  names(sig) <- row.names(dea)
  vi <- viper::msviper(sig, regulon=TS2regulon(as.data.frame(TS)), minsize=minSize, pleiotropy=pleiotropy, verbose=FALSE)
  vi2 <- DataFrame(vi$es[c("nes","size","p.value","nes.bt")])
  colnames(vi2)[3] <- "pvalue"
  vi2$FDR <- p.adjust(vi2$pvalue)
  vi2 <- vi2[order(vi2$pvalue),]
  #vi2$miRNAs <- sapply(row.names(vi2), fam=metadata(TS)$families, FUN=function(x, fam) names(fam)[which(fam==x)])
  vi2
}


#' plMod
#'
#' @param dea A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#' @param var The independent variable to use, either 'sites' or 'score'.
#' @param correctForLength Logical; whether to correct for UTR size (using total bindings as a proxy). 
#' Defaults to TRUE if `var='sites'`, FALSE otherwise.
#'
#' @return a data.frame
#' @export
plMod <- function(dea, TS, minSize=5, var="sites", correctForLength=(var=="sites"), useFC=TRUE){
  print("plMod...")
  dea <- cleanDEA(dea)
  if(useFC){
    sig <- dea$logFC
  } else {
    sig <- -log10(dea$FDR)*sign(dea$logFC)
  }
  names(sig) <- row.names(dea)
  TS <- as.data.frame(TS)
  TS <- aggregate(TS[,c("sites","score")],by=list(family=TS$family, feature=TS$feature),FUN=function(x){ if(is.numeric(x)) return(max(x,na.rm=T)); x[[1]] })
  if(correctForLength){
    ag <- aggregate(TS$sites,by=list(feature=TS$feature),FUN=sum)
    cfl <- ag[,2]
    names(cfl) <- ag[,1]
    cfl <- cfl[names(sig)]
    cfl[which(is.na(cfl))] <- 0
    names(cfl) <- names(sig)
  }else{
    cfl <- NULL
  }
  TS$family <- as.character(TS$family)
  res <- t(sapply(split(TS,TS$family), sig=sig, minSize=minSize, cfl=cfl, FUN=function(x, sig, minSize, cfl){
    r <- c(x$family[1],NA,NA)
    if(nrow(x)<minSize) return(r)
    x <- x[!duplicated(x),]
    row.names(x) <- x$feature
    x2 <- x[names(sig),var]
    x2[which(is.na(x2))] <- 0
    if(is.null(cfl)){
      mod <- try(lm(sig~x2+0),silent=T)
    }else{
      mod <- try(lm(sig~x2+cfl+0),silent=T)
    }
    if(!is(mod,"try-error")) r[2:3] <- c(coef(mod)["x2"], summary(aov(mod))[[1]]["x2","Pr(>F)"])
    return(r)
  }))
  colnames(res) <- c("family","logFC","pvalue")
  res <- DataFrame(res)
  res$logFC <- as.numeric(as.character(res$logFC))
  res$pvalue <- as.numeric(as.character(res$pvalue))
  res$FDR <- p.adjust(res$pvalue)
  row.names(res) <- res$family
  res[order(res$FDR,res$pvalue),-1]
}

#' KS
#'
#' miRNA targets enrichment analysis using a Kolmogorov-Smirnov test on the targets' foldchanges. As all alternatives are considered, significance might not always be consistent with differential miRNA activity.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#'
#' @return a data.frame.
#'
#' @export
KS <- function(DEA, TS, minSize=5, useFC=TRUE){
  print("KS...")
  library(data.table)
  DEA <- cleanDEA(DEA)
  TS$family <- droplevels(TS$family)
  res <- lapply(split(TS,TS$family), DEA=DEA, FUN=function(x,DEA){
    set <- unique(as.character(x$feature))
    set <- intersect(set,row.names(DEA))
    if(useFC){
      targets <- DEA[set,"logFC"]
      nontargets <- DEA[setdiff(row.names(DEA),set),"logFC"]
    } else {
      targets <- -log10(DEA[set,"FDR"]) * sign(DEA[set,"logFC"])
      nontargets <- -log10(DEA[setdiff(rownames(DEA),set),"FDR"]) * sign(DEA[setdiff(row.names(DEA),set),"logFC"])
    }
    ks <- suppressWarnings(try(ks.test(targets, nontargets)$p.value, silent=T))
    if(is(ks,"try-error")) ks <- NA
    list(family=as.character(x$family[1]),
         annotated=length(set),
         ks.pvalue=ks)
  })
  res <- DataFrame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$ks.pvalue,method="fdr")
  row.names(res) <- res$family
  res[order(res$FDR),-1]
}


#' KS2
#'
#' miRNA targets enrichment analysis using a Kolmogorov-Smirnov test on the targets' foldchanges, treating upregulated and downregulated genes separately, and applying a one-sided test to each. Significance here is more likely to incidate an effect consistent with miRNA activity than in a normal KS-test.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#'
#' @return a data.frame.
#'
#' @export
KS2 <- function(DEA, TS, minSize=5, useFC=TRUE){
  print("KS2...")
  library(data.table)
  DEA <- cleanDEA(DEA)
  res <- lapply(split(TS,TS$family), DEA=DEA, FUN=function(x,DEA){
    set <- unique(as.character(x$feature))
    d1 <- DEA[which(DEA$logFC>0),]
    d2 <- DEA[which(DEA$logFC<0),]
    set1 <- intersect(set,row.names(d1))
    set2 <- intersect(set,row.names(d2))
    
    if(useFC){
      targets1 <- d1[set1,"logFC"]
      nontargets1 <- DEA[setdiff(row.names(d1),set1),"logFC"]
      targets2 <- d2[set2,"logFC"]
      nontargets2 <- DEA[setdiff(row.names(d2),set2),"logFC"]
    } else {
      targets1 <- -log10(d1[set1,"FDR"]) * sign(d1[set1,"logFC"])
      nontargets1 <- -log10(DEA[setdiff(rownames(d1),set1),"FDR"]) * sign(DEA[setdiff(row.names(d1),set1),"logFC"])
      targets2 <- -log10(d2[set2,"FDR"]) * sign(d2[set2,"logFC"])
      nontargets2 <- -log10(DEA[setdiff(rownames(d2),set2),"FDR"]) * sign(DEA[setdiff(row.names(d2),set2),"logFC"])
    }
    ks1 <- suppressWarnings(try(ks.test(targets1, nontargets1, alternative="less")$p.value, silent=T))
    ks2 <- suppressWarnings(try(ks.test(targets2, nontargets2, alternative="greater")$p.value, silent=T))
    if(is(ks1,"try-error")) ks1 <- NA
    if(is(ks2,"try-error")) ks2 <- NA
    list(family=as.character(x$family[1]),
         annotated=length(set),
         ks.pvalue.down=ks2,
         ks.pvalue.up=ks1
    )
  })
  res <- DataFrame(rbindlist(res))
  if(!(nrow(res)>1)) return(res)
  res <- res[which(res[,"annotated"]>=minSize),]
  res$ks.pvalue.down <- unlist(res$ks.pvalue.down)
  res$ks.pvalue.up <- unlist(res$ks.pvalue.up)
  res$FDR <- apply(matrix(p.adjust(as.numeric(as.matrix(res[,c("ks.pvalue.down","ks.pvalue.up")])),method="fdr"),ncol=2),1,FUN=min)
  row.names(res) <- res$family
  res[order(res$FDR,apply(res[,grep("pvalue",colnames(res))],1,FUN=min)),-1]
}


#' MW
#'
#' miRNA targets enrichment analysis using a Mann-Whitney / Wilcoxon test on the targets' foldchanges.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of targets for a miRNA family to be tested (default 5).
#'
#' @return a data.frame.
#'
#' @export
MW <- function(DEA, TS, minSize=5, useFC=TRUE){
  print("MW...")
  library(data.table)
  DEA <- cleanDEA(DEA)
  res <-lapply(split(TS,TS$family), DEA=DEA, FUN=function(x,DEA){
    set <- unique(as.character(x$feature))
    set <- intersect(set,row.names(DEA))
    if(useFC){
      targets <- DEA[set,"logFC"]
      nontargets <- DEA[setdiff(row.names(DEA),set),"logFC"]
    } else {
      targets <- -log10(DEA[set,"FDR"]) * sign(DEA[set,"logFC"])
      nontargets <- -log10(DEA[setdiff(rownames(DEA),set),"FDR"]) * sign(DEA[setdiff(row.names(DEA),set),"logFC"])
    }
    mw <- try(wilcox.test(targets, nontargets)$p.value, silent=T)
    if(is(mw,"try-error")) mw <- NA
    list(family=as.character(x$family[1]),
         annotated=length(set),
         wilcox.pvalue=mw)
  })
  res <- DataFrame(rbindlist(res))
  res <- res[which(res[,"annotated"]>=minSize),]
  res$FDR <- p.adjust(res$wilcox.pvalue,method="fdr")
  row.names(res) <- res$family
  res[order(res$FDR),-1]
}


#' gsea
#'
#' miRNA target enrichment analysis among differentially-expressed genes using GSEA.
#'
#' @param DEA A data.frame of the results of a differential expression analysis, with features (e.g. genes) as row names and with at least the following columns: `logFC`, `FDR`
#' @param TS A data.frame of miRNA targets, with at least the following columns: `family`, `rep.miRNA`, `feature`, `sites`.
#' @param minSize The minimum number of elements in a set to be tested (default 5).
#' @param maxSize The maximum number of elements in a set to be tested (default 500). If the test takes too long to run, consider setting this.
#' @param fdr.thres The FDR threshold below which genes are considered; default 0.2.
#' @param nperm The number of permutations, default 2000. The more permutations, the more precise the p-values, in particular for small p-values.
#'
#' @return a data.frame.
#'
#' @importFrom fgsea fgsea
#' @export
gsea <- function(DEA, TS, minSize=5, maxSize=300, fdr.thres=0.5, nperm=2000, useFC=TRUE){
  print("GSEA...")
  DEA <- cleanDEA(DEA)
  DEA <- DEA[which(DEA$FDR<=fdr.thres),]
  if(useFC){
    sig <- DEA$logFC
  } else {
    sig <- -log10(DEA$FDR)*sign(DEA$logFC)
  }
  names(sig) <- row.names(DEA)
  sets <- lapply(split(TS$feature,TS$family),tested=names(sig),FUN=function(x,tested){ intersect(unique(x),tested) })
  sets <- sets[which(sapply(sets,length)>=minSize)]
  res <- fgsea::fgsea(sets, sig, nperm, minSize=minSize, maxSize=maxSize)
  res <- DataFrame(res[order(res$padj,res$pval),])
  colnames(res)[1:5] <- c("family","pvalue","FDR","ES","normalizedEnrichment")
  colnames(res)[8] <- "features"
  row.names(res) <- res$family
  return(res[,-1])
}
