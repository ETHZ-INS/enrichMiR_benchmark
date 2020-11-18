# Data  
- **Summary**  
Inputs for the `tgermade_thesis` & other analyses.  
- **Content**   
DEAs (SummarizedExperiment), TargetScan annotations (DataFrame) & other relevant files  
- **Input for**   
`../functions/enrichMiR_doRun.R`; miscellaneous  
- **Generated via**  
`DEA.SE.rds` files: `../*datagen.DEA.SE.Rmd`; TargetScan files: `../functions/getTS.R`; miscellaneous  
- **Additional info**  
Contains DEA files for benchmark datasets amin, bartel, cherone, jeong, ratPolyA, whipple.  
Contains TargetScan files filtered for species human, mouse, rat.  
Contains set of genes that were excluded from bartel DEAs (`bartel.hek.excludedGenes.rds` & `bartel.hela.excludedGenes.rds`)
