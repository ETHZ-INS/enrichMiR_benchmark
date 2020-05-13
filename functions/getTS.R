#' getTS
#'
#' @param species character object. Can be "human", "mouse" or "rat"
#'
#' @return TargetScan miRNA target dataframe with family information in metadata()
#'
getTS <- function(species=c("human","mouse","rat")){
  library(S4Vectors)
  species <- match.arg(species)
  
  # assign species ID
  spec <- switch( species,
                  human = 9606,
                  mouse = 10090,
                  rat = 10116 )
  
  
  # download TargetScan miRNA targeting dataset
  tmp <- tempfile()
  if(species=="human"){
    download.file(
      "http://www.targetscan.org/vert_72/vert_72_data_download/Summary_Counts.all_predictions.txt.zip",tmp)
    unzip(tmp, exdir="data/")
    full <- fread("data/Summary_Counts.all_predictions.txt")
  } else if(any(species %in% c("mouse","rat"))){
    download.file(
      "http://www.targetscan.org/mmu_72/mmu_72_data_download/Summary_Counts.all_predictions.txt.zip", tmp)
    unzip(tmp, exdir="data/")
    full <- fread("data/Summary_Counts.txt")
  }
  
  # limit to selected species
  sub <- full[full$'Species ID' == spec,]
  
  # generate TargetScan dataframe
  ts <- data.frame(family = sub$'miRNA family',
                   rep.miRNA = sub$'Representative miRNA',
                   feature = sub$'Gene Symbol',
                   sites = sub$'Total num conserved sites',
                   score = as.numeric(as.character(sub$'Cumulative weighted context++ score'))
  )
  # aggregate
  ts <- DataFrame(
    aggregate(ts[,c("sites","score")], by=ts[,c("family","feature","rep.miRNA")], na.rm=TRUE, FUN=mean)
  )
  
  # download TargetScan miRNA families dataset
  tmp <- tempfile()
  if(species=="human"){
    download.file(
      "http://www.targetscan.org/vert_72/vert_72_data_download/miR_Family_Info.txt.zip",tmp)
  } else if(any(species %in% c("mouse","rat"))){
    download.file(
      "www.targetscan.org/mmu_72/mmu_72_data_download/miR_Family_Info.txt.zip", tmp)
  }
  unzip(tmp, exdir="data/")
  full <- fread("data/miR_Family_Info.txt")
  
  # limit to selected species
  sub <- full[full$'Species ID' == spec,]
  
  fam <- sub$`Seed+m8`
  names(fam) <- sub$`MiRBase ID`
  
  # add family info to ts dataframe as attribute
  metadata(ts)$families <- fam
  
  # enrichMiR cant handle 0 values for sites feature
  return(ts[ts$sites!=0,])
}