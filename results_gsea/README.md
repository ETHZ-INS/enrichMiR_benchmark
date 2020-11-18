# GSEA enrichMiR test results  
- **Summary**  
EnrichMiR test results for comparing performances when using different inputs & primary tests.
- **Content**  
lists of dataframes of datasets amin, bartel, cherone, jeong, ratPolyA
- **Analyzed in**  
`../analysis/gsea_variants.Rmd` (html available)
- **Generated via**  
`../functions/gsea_variants_doRun.R`
- **Additional info**  
  - For each run different available & possible inputs were supplied to the tests (*signal*). The tests were either run in their original / pure form (w/out primary test) or with supplied results from other enrichMiR tests (w/ primary test). The supplied tests stem from `../results/`.    
  - Variables (miRNAs) for the secondary test were defined by selecting the top 10% of instances in the primary test.  
  - The naming system is the following: `dataset.gsea.signaltype.rds` & `dataset.gsea.signaltype_feature.rds` in cases where `signaltype` is continuous:   
    `signaltype`: bin (binary) / cont (continuous) ;  
    `feature`: fc (logFC signal) / fdr (signed_logFDR signal)     

The files are stored on a local server under `/mnt/schratt/enrichMiR_benchmark/results_gsea/`.
