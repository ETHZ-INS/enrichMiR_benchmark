completeEnrichMiR <- function(p){
  
  # enrichMiR
  
  source("functions/runEnrichMiR.R")
  
  ## extract DEAs from SE rowData and load annotation
  p <- loadAll(p)
  ## clean up DEA lists
  p <- cleanDEA(p)
  ## generate annotation library permutations based on props
  if(p$perm)  p <- TSPerm(p)
  ## run enrichMiR
  # e.list <- getEnrichMiR(p, cores=1)
  e.list <- readRDS("e.list.rds")
  ## aggregate some test results
  test.list <- combineTests(e.list, tests=NULL, fisher.agg=TRUE, geom.agg=TRUE)

  # benchmark
  
  source("functions/runBenchmark.R")
  
  ## get benchmark results dataframe
  bm.df <- getBench(test.list, p$tp, p$perm)
  
  return(bm.df)
}