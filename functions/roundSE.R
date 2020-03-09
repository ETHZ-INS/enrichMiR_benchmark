######### reduce dataset size function ######### 
# tgermade
# 05.02.20


#' roundSE
#'
#' @param se SummarizedExperiment object
#'
#' @return SE object with rounded digits
#'
roundSE <- function(se){
  
  library(SummarizedExperiment)
  
  # round every assay except for 'counts'
  for(a in setdiff(assayNames(se),"counts")){
    for(i in 1:ncol(se)){ 
      assays(se)[[a]][,i] <- plgINS::dround(assays(se)[[a]][,i], 
                                            roundGreaterThan1 = TRUE)
    }
    assays(se)[[a]][is.na(assays(se)[[a]])] <- 0
  }
  # round every rowData column
  for(i in 1:ncol(rowData(se))){
    col <- rowData(se)[,i]
    ## round dataframe columns
    if(is.data.frame(col)){
      for(j in 1:ncol(col)){
        if(is.numeric(col[[j]])){
          col[[j]] <- plgINS::dround(col[[j]], roundGreaterThan1 = TRUE)
          col[[j]][is.na(col[[j]])] <- 0
        }
      }
    }
    ## round normal columns
    if(is.numeric(col)){
      col <- plgINS::dround(col, roundGreaterThan1 = TRUE)
      col[is.na(col)] <- 0
    }
    rowData(se)[[i]] <- col
  }
  return(se)
}
