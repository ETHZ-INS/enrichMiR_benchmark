# source("../master_19_20/enrichMiR/enrichMiR/R/tests.R")
# 
# e.list <- readRDS("results2/bartel.hek.enrichMiR.gene.rds")
# e <- e.list$miR.122@res$modscore
# 
# th.abs.logFC=0
# th.FDR=0.05
# binary.signatures <- list(
#   down=.dea2binary(x, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=-1),
#   up=.dea2binary(x, th=th.FDR, th.alfc=th.abs.logFC, restrictSign=1)
# )
# signal <- binary.signatures$down
# bm <- sapply(split(as.character(sets$feature), sets$set), FUN=function(x) names(signal) %in% x)
# bm <- bm[,colSums(bm)>0]
# 
# co <- rownames(e[e$FDR<=.05,])
# 
# signald <- data.frame( y=signal, bm[,co,drop=FALSE] )
# 
# form <- y~0+.
# mod <- glm.cons( form, data=signald, family="binomial", cons=1, cons.inter=-1)
# 
# 
# res <- coef(summary(mod))
# res <- res[order(res[,4]),,drop=FALSE]
# colnames(res) <- c("beta","stderr","z","pvalue")
# res <- res[grep("^\\(Intercept\\)$|FALSE$", row.names(res), invert=TRUE),,drop=FALSE]
# row.names(res) <- gsub("TRUE","",row.names(res))
# 
# head(res)
# 
# tp <- c(miR.122 = "GGAGUGU")
# which(rownames(res)==tp)
# nrow(res)
# 
# 
# ########################################################################
# 
# 
# 
# data <- signald
# data$y <- rnorm(nrow(data))
# colnames(data) <- colnames(signald)
# 
# test1 <- glm(form, data, family="gaussian")
# test2 <- lm(form, data)



########################################################################

devtools::load_all("../master_19_20/enrichMiR/enrichMiR/")

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


#' regmir.new
#'
#' @param e 
#' @param sets 
#' @param dea 
#' @param binary 
#' @param alpha 
#' @param use.intercept 
#' @param keepAll 
#' @param cons 
#' @param KO 
#'
#' @return
#' @export
#'
#' @examples
regmir.new <- function(e, sets, dea, sig.binary=FALSE, var.binary=FALSE, cpm.weights=TRUE,
                       use.intercept=FALSE, keepAll=TRUE, coeff.cons=TRUE, KO=TRUE, useFC=TRUE){
  
  suppressPackageStartupMessages(c(
    library(glmnet),
    library(zetadiv)
  ))

  e <- as.data.frame(e)
  e[is.na(e)] <- 1
  
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
      signal <- -log10(dea$FDR)*sign(dea$logFC)
    }
    names(signal) <- rownames(dea)
  }
  
  # gnerate bm
  ## prepare the target matrix
  if(var.binary){
    bm <- sapply(split(as.character(sets$feature), sets$set), FUN=function(x) names(signal) %in% x)
  }else{
    TS2 <- aggregate(sets[,"score",drop=FALSE], by=as.data.frame(sets[,c("set","feature")]), FUN=min)
    TS2 <- split(TS2[,c("score","feature")], TS2$set)
    bm <- sapply(TS2, FUN=function(x){
      row.names(x) <- x$feature
      -1*x[names(signal),"score"]
    })
    bm[is.na(bm)] <- 0
    colnames(bm) <- names(TS2)
  }
  bm <- bm[,colSums(bm)>0]
  
  # generate co
  if(sum(e$FDR<=.05)>0){
    co <- rownames(e[e$FDR<=.05,])
  } else if(sum(e$FDR<=.2)>0){
    co <- rownames(e[e$FDR<=.2,])
  } else {
    co <- rownames(e[1:10,])
  }
  # combine signal, bm and co
  signald <- data.frame( y=signal, bm[,co,drop=FALSE] )
  
  # new fit to get significance estimates
  if(use.intercept){
    form <- y~.
  }else{
    form <- y~0+.
  }
  if(cpm.weights){
    w <- dea[rownames(signald),]$logCPM + abs(min(dea[rownames(signald),]$logCPM))
    if(sig.binary){
      if(coeff.cons){
        mod <- glm.cons( form, data=signald, family="binomial", cons=1, cons.inter=-1, weights=w)
      } else {
        mod <- glm( form, data=signald, family="binomial", weights=w )
      }
    }else{
      mod <- lm( form, data=signald, weights=w )
    }
  } else {
    if(sig.binary){
      if(coeff.cons){
        mod <- glm.cons( form, data=signald, family="binomial", cons=1, cons.inter=-1)
      } else {
        mod <- glm( form, data=signald, family="binomial" )
      }
    }else{
      mod <- lm( form, data=signald )
    }
  }
  
  # we extract the coefficients and p-values, and reorganize the output:
  res <- coef(summary(mod))
  res <- res[order(res[,4]),,drop=FALSE]
  colnames(res) <- c("beta","stderr",ifelse(sig.binary,"z","t"),"pvalue")
  res <- res[grep("^\\(Intercept\\)$|FALSE$", row.names(res), invert=TRUE),,drop=FALSE]
  row.names(res) <- gsub("TRUE","",row.names(res))
  
  if(keepAll){
    co2 <- setdiff(colnames(bm), co)
    co2 <- data.frame(row.names=co2, beta=rep(NA_real_,length(co2)),
                      stderr=NA_real_, z=NA_real_, pvalue=1)
    if(!sig.binary) colnames(co2)[3] <- "t"
    res <- rbind(res,co2)
  }
  
  # we adjust using all features as number of comparisons
  if(nrow(res)>0){
    res$FDR <- p.adjust(res$pvalue, n=ncol(bm))
  }
  res
}

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
dea <- loadAll("../master_19_20/enrichMiR_benchmark/data/bartel.hek.DEA.SE.rds",TPs)
dea <- lapply(dea, cleanDEA)
e.list <- readRDS("../master_19_20/enrichMiR_benchmark/results/bartel.hek.enrichMiR.rds")
TS <- readRDS("../master_19_20/enrichMiR_benchmark/data/TargetScan_human.rds")
if(all(c("family", "feature") %in% colnames(TS))){
  colnames(TS) <- gsub("^family$","set", colnames(TS))
}
TS <- DataFrame(TS)




for(w in cpm.weights <- c(T,F)){
  print(paste("weights:", w))
  W=ifelse(w,"w","nw")
  
  for(v in var.bin <- c(T,F)){
    print(paste("binary variables:", v))
    V=ifelse(v,"bin","cont")
    
    for(s in sig.bin <- c(T,F)){
      print(paste("binary signal:", s))
      S=ifelse(s,"bin","cont")
      
      if(s){
        for(c in coeff.cons <- c(T,F)){
          print(paste("constrained coefficients:", c))
          C=ifelse(c,"c","nc")
          
          e <- lapply(names(dea), function(i) lapply(e.list[[i]]$original@res, function(e) 
            regmir.new(e, TS, dea[[i]], KO=FALSE, cpm.weights=w, var.binary=v, sig.binary=s, coeff.cons=c) 
            ))
          saveRDS(e, paste0("results2/bartel.hek.regmir.", S,"_",V,".",W,".",C,".rds") )
          print( paste(paste0("bartel.hek.regmir.", S,"_",V,".",W,".",C,".rds"),"saved.") )
          
        }
      } else {

        e <- lapply(names(dea), function(i) lapply(e.list[[i]]$original@res, function(e) 
          regmir.new(e, TS, dea[[i]], KO=FALSE, cpm.weights=w, var.binary=v, sig.binary=s) 
          ))
        saveRDS(e, paste0("results2/bartel.hek.regmir.", S,"_",V,".",W,".rds") )
        print( paste(paste0("bartel.hek.regmir.", S,"_",V,".",W,".rds"),"saved.") )
        
      }
    }
  }
}




# for(w in cpm.weights <- c(T,F)){
#   print(paste("weights:", w))
#   W=ifelse(w,"w","nw")
#   
#   for(v in var.bin <- c(T,F)){
#     print(paste("binary variables:", v))
#     V=ifelse(v,"bin","cont")
#     
#     for(s in sig.bin <- c(T,F)){
#       print(paste("binary signal:", s))
#       S=ifelse(s,"bin","cont")
#       
#       if(s){
#         for(c in coeff.cons <- c(T,F)){
#           print(paste("constrained coefficients:", c))
#           C=ifelse(c,"c","nc")
#           
#           
#           print(paste0("results/amin.regmir.", S,"_",V,".",W,".",C,".rds") )
#           print( paste(paste0("results/amin.regmir.", S,"_",V,".",W,".",C,".rds"),"saved.") )
#         }
#       } else {
# 
#         print(paste0("results/amin.regmir.", S,"_",V,".",W,".rds") )
#         print( paste(paste0("results/amin.regmir.", S,"_",V,".",W,".rds"),"saved.") )
#       }
#     }
#   }
# }
# 
# 
# 
# 
# 
# 
# 
# 
# 
# 
# # do logistic GLM with coeff constr
# print("logistic run...")
# e.logistic <- lapply(names(dea), function(i) lapply(e.list[[i]]$original@res, function(e){
#   print("testing...")
#   regmir.new(e, TS, dea[[i]], KO=TRUE, binary=TRUE) }))
# saveRDS(e.logistic, "results2/regmir.amin.log.rds")
# 
# # do logistic GLM without coeff constr
# print("logistic no-constraint run...")
# e.logistic.nc <- lapply(names(dea), function(i) lapply(e.list[[i]]$original@res, function(e) 
#   regmir.new(e, TS, dea[[i]], KO=TRUE, binary=TRUE, cons=FALSE) ))
# saveRDS(e.logistic.nc, "results2/regmir.amin.lognc.rds")
# 
# # do non-logistic GLM
# print("non-logistic run...")
# e.gaussian <- lapply(names(dea), function(i) lapply(e.list[[i]]$original@res, function(e) 
#   regmir.new(e, TS, dea[[i]], KO=TRUE, binary=FALSE) ))
# saveRDS(e.gaussian, "results2/regmir.amin.nonlog.rds")
# 
# 
# 
# log <- readRDS("results2/regmir.amin.log.rds")
# lognc <- readRDS("results2/regmir.amin.lognc.rds")
# nonlog <- readRDS("results2/regmir.amin.nonlog.rds")
# 
# names(log) <- names(lognc) <- names(nonlog) <- names(e.list)
# 
# d1 <- data.frame(p.rank=sapply(e.list$`218DKOvWT`$original@res, function(x) which(rownames(x)==TPs[names(log)])),
#            p.sig=sapply(e.list$`218DKOvWT`$original@res, function(x) x$FDR[rownames(x)==TPs[names(log)]]<=.05),
#            orig="original")
# d2 <- data.frame(p.rank=sapply(log$`218DKOvWT`, function(x) which(rownames(x)==TPs[names(log)])),
#            p.sig=sapply(log$`218DKOvWT`, function(x) x$FDR[rownames(x)==TPs[names(log)]]<=.05),
#            orig="regmir.logistic")
# d3 <- data.frame(p.rank=sapply(lognc$`218DKOvWT`, function(x) which(rownames(x)==TPs[names(lognc)])),
#            p.sig=sapply(lognc$`218DKOvWT`, function(x) x$FDR[rownames(x)==TPs[names(log)]]<=.05),
#            orig="regmir.logistic.nc")
# d4 <- data.frame(p.rank=sapply(nonlog$`218DKOvWT`, function(x) which(rownames(x)==TPs[names(nonlog)])),
#            p.sig=sapply(nonlog$`218DKOvWT`, function(x) x$FDR[rownames(x)==TPs[names(log)]]<=.05),
#            orig="regmir.nonlogistic")
# d <-rbind(d1,d2,d3,d4)
# d$method <- rownames(d)[1:nrow(d1)]
# for(x in colnames(d)[-1]){
#   d[[x]] <- as.factor(d[[x]])
# }
# d$method <- droplevels(d$method)
# 
# ggplot(d,aes(method,p.rank,color=orig,shape=p.sig)) + geom_point(size=3) + scale_y_sqrt() + coord_flip()
# 
# 
# d1 <- getBench(e.list, TPs[names(log)], perm = TRUE)
# d2 <- getBench(log, TPs[names(log)], e.obj = FALSE)
# d3 <- getBench(lognc, TPs[names(log)], e.obj = FALSE)
# d4 <- getBench(nonlog, TPs[names(log)], e.obj = FALSE)
# d1 <- d1[d1$prop.rep=="original",-c(2,9)]
# 
# df <- dplyr::bind_rows(list(original=d1,logistic=d2,logistic.nc=d3,nonlogistic=d4), .id="orig")
# 
# TPvFP(df, label.all=TRUE, features = "orig")
# HM(df, form=method~orig)
# HM(df, metric="FP.atFDR05", form=method~orig)
