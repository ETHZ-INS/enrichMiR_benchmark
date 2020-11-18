# EnrichMiR test results  
- **Summary**  
EnrichMiR test results using the enrichMiR branch freeze20200803.  
- **Content**  
`enrichMiR.rds` (enrich.results) & `benchmark.rds` (dataframe) files for datasets amin, bartel, cherone, jeong, ratPolyA  
- **Analyzed in**  
`../tgermade_thesis/tgermade_thesis_Rcode.Rmd`  
- **Generated via**  
`../functions/enrichMiR_doRun.R`  
- **Additional info**  
  - For each dataset there are results of running enrichMiR 1) on `DEA.SE.rds` files & w/out mirexpr (species-specific miRNA background), 2) on `exon.DEA.SE.rds` files & w/out mirexpr, 3) on `DEA.SE.rds` files & w/ mirexpr (tissue-specific miRNA background). The DEA files can be found under `../data/`.    
  - The naming system is the following: `dataset.mode.benchmark.rds` & `dataset.mode.enrichMiR.rds`:   
    `mode`: - (empty for standard run) / exon (exon-specific DEA) / mirexpr (tissue-specific background) 
