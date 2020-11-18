# Continuous enrichMiR test results  
- **Summary**  
EnrichMiR test results for comparing performances when using 2 different types of signal input (either logFC or signed_logFDR).  
- **Content**  
enrichMiR.rds (enrich.results) & benchmark.rds (dataframe) files for datasets amin, bartel, cherone, jeong, ratPolyA  
- **Analyzed in**  
`../analysis/enrichMiR_results.continuous.Rmd` (html available)  
- **Generated via**  
`../functions/cont_variants_doRun.R`  
- **Additional info**  
The signed_logFDR feature is defined as `gn(logFC)*-log10(FDR)`.  

The files are stored on a local server under `/mnt/schratt/enrichMiR_benchmark/results_continuous/`.
