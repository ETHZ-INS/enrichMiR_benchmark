######### Differential Expression Analysis function #########
# tgermade & cschlaeppi
# 27.04.20

#' DEA
#'
#' @param se SummarizedExperiment object
#' @param use An id vector to subset 'se': either c(TRUE,FALSE,...) or c(1:9). Default: no subsetting
#' @param name Suffix used for DEA results in rowData
#' @param model A formula for linear regression. Model to be compared to 'model0'.
#' @param model0 A formula for linear regression. Model to compare 'model' to. Default: ~1
#' @param testCoeff Select which model ids (columns) to use for F-test. Default: difference btw. 'model' & 'model0'
#' @param QLF Logical value if QLF should be applied or not (default = FALSE). Don't use QLF if low replicate number.
#'
#' @return A SummarizedExperiment with additional DEA results dataframe in rowData(), additional 
#'          logFC assay and logcpm assay (if not already present)
#'
#' @examples se <- DEA(se, use = se$stage %in% c(0,9), name = "0v9", model = ~stage, 
#'                  testCoeff = c("stage0","stage9"))
#'
DEA <- function(se, use = NULL, name, model, model0 = ~1, testCoeff = NULL, QLF = FALSE){
  
  library(edgeR)
  
  # allocation
  
  ## se subset & drop unused factor levels
  if( !is.null(use) ){
    se.sub <- se[,use]
  } else { se.sub <- se }
  
  colData(se.sub) <- droplevels(colData(se.sub))
  
  ## models
  mm <- model.matrix(model, data=as.data.frame(colData(se.sub)))
  mm0 <- model.matrix(model0, data=as.data.frame(colData(se.sub)))
  
  if( is.null(testCoeff) ){
    testCoeff <- setdiff(colnames(mm), colnames(mm0))
  }
  
  # DEA  
  
  ## create 'DGEList' object
  comb_dds <- DGEList(assay(se.sub))
  ## To do the normalization
  comb_dds <- calcNormFactors(comb_dds)
  ## estimate dispersion
  comb_dds <- estimateDisp(comb_dds, mm)
  
  if (QLF) {
    ## fit negative binomial distribution on counts per gene
    comb_fit <- glmQLFit(comb_dds, mm)
    ## fit a GLM on the data
    comb_lrt <- glmQLFTest(comb_fit, testCoeff)
  } else {
    ## fit negative binomial distribution on counts per gene
    comb_fit <- glmFit(comb_dds, mm)
    ## fit a GLM on the data
    comb_lrt <- glmLRT(comb_fit, testCoeff)
  }
  
  ## top genes that change relative to stage 0
  comb_res <- as.data.frame(topTags(comb_lrt, Inf))
  
  ## add symbol names to DEA results
  comb_res <- comb_res[rownames(se.sub),]
  # comb_res$symbol <- rowData(se)$symbol
  
  # add DEA results to SE
  dea_name <- paste0("DEA.", name)
  rowData(se)[[dea_name]] <- DataFrame(comb_res)
  
  return(se)
}
