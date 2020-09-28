# Goal: for tests that use continuous signal input, should we use logFC or signed logFDR?
# enrichMiR: uses branch freeze20200803 with enrichMiR.sub.R & tests.sub.R
# output: enrichMiR results of tests using continuous signals for all benchmark datasets 

print("loading & setup")

suppressPackageStartupMessages({
  library(SummarizedExperiment)
  library(SEtools)
  library(plgINS)
  library(edgeR)
  library(DT)
  library(viridis)
  library(data.table)
  library(dplyr)
  library(tibble)
  library(matrixStats)
  library(RUVSeq)
  library(BiocParallel)
})
devtools::load_all("/mnt/schratt/tgermade_test/hiwi_20/enrichMiR.vfreeze2/enrichMiR/")


# miRNAs (mirexpr vector generation)
getMirexpr <- function(tissue){
  library(microRNAome)
  data("microRNAome")
  counts <- rowMedians(cpm(calcNormFactors(DGEList(assay(microRNAome)[,microRNAome$cell_tissue==tissue]))))
  names(counts) <- rownames(microRNAome)
  return(counts[counts>=10])
}

## human
hela <- getMirexpr("hela")
hela <- c(hela, structure(rep(1000,16), names=c(
  "hsa-let-7a-5p","hsa-miR-1-3p","hsa-miR-124-3p.1","hsa-miR-137","hsa-miR-139-5p",
  "hsa-miR-143-3p","hsa-miR-144-3p","hsa-miR-153-3p","hsa-miR-155-5p","hsa-miR-182-5p",
  "hsa-miR-199a-5p","hsa-miR-204-5p","hsa-miR-205-5p","hsa-miR-216b-5p","hsa-miR-223-3p",
  "hsa-miR-7-5p")) )
hek <- getMirexpr("hek293")
hek <- c(hek, structure(rep(1000,12), names=c(
  "hsa-miR-122-5p","hsa-miR-133a-3p.1","hsa-miR-138-5p","hsa-miR-145-5p",
  "hsa-miR-184","hsa-miR-190a-5p","hsa-miR-200b-3p","hsa-miR-216a-5p",
  "hsa-miR-217","hsa-miR-219a-5p","hsa-miR-375","hsa-miR-451a")) )
## rat
ratmirna <- read.csv("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/DIV20_Neuron_HC_miRNA.csv")
rat <- 2^(ratmirna$logCPM) 
names(rat) <- ratmirna$Gene.name
rat <- rat[rat>=10]
## mouse
mousemirna <- readRDS("/mnt/schratt/tgermade_test/master_19_20/enrichMiR_benchmark/data/miRNA_celltype_mouseBrain_GSE30286.SE.rds")
mouse <- rowMedians(cpm(calcNormFactors(DGEList(assay(mousemirna)[,mousemirna$type!="Cerebellum"]))))
names(mouse) <- rownames(mousemirna)
mouse <- mouse[mouse>=10]
mouse.df <- data.frame(expression=mouse, id=sub("\\*","",names(mouse)), row.names=names(mouse))
mouse.df <- aggregate(mouse.df$expression, by=list(mouse.df$id), FUN="max")
mouse.df <- data.frame(expression=rep(mouse.df$x,3),
                       id=c(as.character(mouse.df$Group.1), 
                            sapply(mouse.df$Group.1, function(x) paste0(x,"-3p")),
                            sapply(mouse.df$Group.1, function(x) paste0(x,"-5p"))) )
mouse <- structure(mouse.df$expression, names=as.character(mouse.df$id))


mirexpr <- list(hela=hela, 
                hek=hek, 
                rat=rat,
                mouse=mouse)


alldea=c("bartel.hela.DEA.SE.rds",
        "bartel.hek.DEA.SE.rds",
        "bartel.hela.exon.DEA.SE.rds",
        "bartel.hek.exon.DEA.SE.rds",
        "jeong.DEA.SE.rds",
        "jeong.exon.DEA.SE.rds",
        "ratPolyA.DEA.SE.rds",
        "ratPolyA.exon.DEA.SE.rds",
        "amin.DEA.SE.rds",
        "amin.exon.DEA.SE.rds",
        "cherone.d0.DEA.SE.rds",
        "cherone.d1.DEA.SE.rds",
        "cherone.d0.exon.DEA.SE.rds",
        "cherone.d1.exon.DEA.SE.rds")

allann <- c(hsa="TargetScan_human.rds",
            mmu="TargetScan_mouse.rds",
            rno="TargetScan_rat.rds")

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

p <- list(
  dea=NULL,
  ann=NULL,
  tp=TPs,
  tests=NULL,
  mirnas=NULL,
  perm=FALSE,
  props=NULL
)

for(dea in alldea){
  print(paste("running enrichMiR on:", dea, "..."))
  # setup p
  p$dea <- paste0("../master_19_20/enrichMiR_benchmark/data/", dea)
  if(grepl("bartel", dea) | grepl("jeong", dea)){
    p$ann <- paste0("../master_19_20/enrichMiR_benchmark/data/", allann[["hsa"]])
  } else if(grepl("rat", dea)){
    p$ann <- paste0("../master_19_20/enrichMiR_benchmark/data/", allann[["rno"]])
  } else {
    p$ann <- paste0("../master_19_20/enrichMiR_benchmark/data/", allann[["mmu"]])
  }
  
  for(mirlib in 0:1){
    if(mirlib){
      if(grepl("bartel.hela", dea)){
        p$mirnas <- mirexpr$hela
      } else if(grepl("bartel.hek", dea)){
        p$mirnas <- mirexpr$hek
      } else if(grepl("rat", dea)){
        p$mirnas <- mirexpr$rat
      } else if(grepl("cherone", dea) | grepl("amin", dea)){
        p$mirnas <- mirexpr$mouse
      }
    }
    
    source("functions/runEnrichMiR.R")
    
    # extract DEAs from SE rowData and load annotation
    p <- loadAll(p)
    # clean up DEA lists
    p <- cleanDEA(p)
    # run enrichMiR
    e.list <- getEnrichMiR(p, cores=1)
    
    if(!mirlib){
      saveRDS(e.list, file=file.path("results", gsub("DEA\\.SE","enrichMiR", dea)))
    } else {
      saveRDS(e.list, file=file.path("results", gsub("DEA\\.SE","mirexpr.enrichMiR", dea)))
    }
    print("successfully processed & saved.")
  }
    
}

