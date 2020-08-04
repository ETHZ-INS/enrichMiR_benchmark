#--------------------------------------------------------------------------------
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


#--------------------------------------------------------------------------------
#' getBench
#'
#' @param e.list list of dataframes containing enrichMiR results
#' @param TP vector of treatment miRNAs seeds and their names
#' @param perm binary; TRUE if enrichMiR was done with annotation permutations
#'
#' @return a dataframe of pooled enrichMiR test scores
#' 
getBench <- function(e.list, TP, perm){ 
  # generate the benchmarking scores
  if(!perm){
    bm.list <- lapply(names(TP), function(x) doBenchmark(e.list[[x]]@res, TP[x]))
    names(bm.list) <- names(TP)
  } else {
    bm.list <- lapply(names(TP), function(x) lapply(e.list[[x]], function(e)
      doBenchmark(e@res, TP[x])))
    for(x in names(TP)){
      names(bm.list[[x]]) <- names(e.list[[x]])
    }
    bm.list <- lapply(bm.list, FUN=function(x) dplyr::bind_rows(x, .id = "prop.rep"))
  }
  # generate a results df for plotting
  bm.df <- dplyr::bind_rows(bm.list, .id="treatment")
  if(perm){
    bm.df$prop <- unlist(lapply(strsplit(bm.df$prop.rep, "[.]"), FUN=function(x) x[1]))
    bm.df$prop <- factor(bm.df$prop, levels=unique(bm.df$prop))
  }
  return(bm.df)
}
