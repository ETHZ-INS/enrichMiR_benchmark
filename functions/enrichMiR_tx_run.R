

#' getEnriched
#'
#' @param deas 
#' @param ann 
#' @param species 
#' @param n.out 
#'
#' @return
#' @export
#'
#' @examples
getEnriched <- function(deas, ann=NULL, species=c("hsa","mmu","rno"), n.out=NULL, tests=NULL){
  devtools::load_all("/mnt/schratt/tgermade_test/master_19_20/enrichMiR/enrichMiR/")
  if(is.null(n.out)){
    if(is.character(deas)){
      n.out <- rev(strsplit(deas, "\\/")[[1]])[1]
    } else n.out <- "enrich.out.rds"
  }
  if(is.character(deas)){
    deas <- readRDS(deas)
  }
  if(is.null(ann)){
    ann <- switch( species,
                    hsa = "/mnt/schratt/miRNA_KD/EnrichMir_Work/Aggregate_fg/hsa.utrs.agg_bartel.fst",
                    mmu = "/mnt/schratt/miRNA_KD/EnrichMir_Work/Aggregate_fg/mmu.utrs.agg_bartel.fst",
                    rno = "/mnt/schratt/miRNA_KD/EnrichMir_Work/Aggregate_fg/rno.utrs.agg_bartel.fst" )
  }
  if(is.character(ann)){
    source("/mnt/schratt/tgermade_test/hiwi_20/scanMiR/scanMiR/R/IndexedFst.R")
    library(fst)
    ann <- loadIndexedFst(ann)
  }
  
  # prep: dea
  deas <- lapply(deas, function(dea){
    rownames(dea) <- sapply(strsplit(row.names(dea), "\\."), FUN=function(x) x[1]) 
    names(dea)[2] <- "logFC"
    names(dea)[4] <- "PValue"
    names(dea)[5] <- "FDR"
    # there can't be any NA FDRs for some of the tests
    dea$FDR[is.na(dea$FDR)] <- 1
    dea
  })
  # prep: annotation
  ann <- as.data.table(as.data.frame(ann))
  ann <- ann[ann$n.8mer>0 | ann$n.7merm8>0 | ann$n.7merA1>0,]
  ann$sites <- ann$n.8mer + ann$n.7merm8 + ann$n.7merA1
  
  # enrichMiR
  e.list <- lapply(deas[1:3], function(dea){
    enrichMiR(DEA=dea, TS=as.data.frame(ann), tests=tests, field="logFC")
    })
  names(e.list) <- names(deas)
  saveRDS(e.list, paste0("results4/",n.out))
}


deas <- "/mnt/schratt/enrichMir_datasets/bartel.salmon/DEAs.shrinked.rds"
getEnriched(deas, species="hsa", n.out="bartel.enrichMiR.rds", tests="regmir")




