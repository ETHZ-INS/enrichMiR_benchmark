# Regmir enrichMiR test results  
- **Summary**  
EnrichMiR test results for comparing performances when using different inputs & primary tests.
- **Content**  
lists of dataframes of datasets amin, bartel, cherone, jeong, ratPolyA
- **Analyzed in**  
`../analysis/regmir_variants.Rmd` (html available)
- **Generated via**  
`../functions/regmir_variants_doRun.R`
- **Additional info**  
  - For each run different available & possible inputs were supplied to the tests (*signal*, *variables*, *weights*, additional options). The tests were either run in their original / pure form (w/ lasso regression as primary test: cv.glmnet) or with supplied results from other enrichMiR tests. The supplied tests stem from `../results/`.    
  - Variables (miRNAs) for the secondary test were defined by selecting significant instances in the primary test.  
  - The naming system is the following: `dataset.regmir.signaltype_variabletype.weights.rds` & `dataset.regmir.signaltype_variabletype.weights.constraint.rds` in cases where signaltype is binary:   
    `signaltype`: bin (binary) / cont (continuous) ;  
    `variabletype`: bin / cont ;  
    `weights`: w (w/ logCPM weights) / nw (w/out logCPM weights) ;  
    `constraint`: c (w/ constrained coefficients) / nc (w/out constrained coefficients)  

The files are stored on a local server under `/mnt/schratt/enrichMiR_benchmark/results_regmir/`.
