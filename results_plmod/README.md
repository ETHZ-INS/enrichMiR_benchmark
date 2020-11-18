# Plmod enrichMiR test results  
- **Summary**  
EnrichMiR test results for comparing performances when using different inputs & primary tests.
- **Content**  
lists of dataframes of datasets amin, bartel, cherone, jeong, ratPolyA
- **Analyzed in**  
`../analysis/plmod_variants.Rmd` (html available)
- **Generated via**  
`../functions/plmod_variants_doRun.R`
- **Additional info**  
  - For each run different available & possible inputs were supplied to the tests (*signal*, *variables*, *weights*). The tests were either run in their original / pure form (w/out primary test) or with supplied results from other enrichMiR tests (w/ primary test). The supplied tests stem from `../results/`.    
  - Variables (miRNAs) for the secondary test were defined by selecting significant instances in the primary test.   
  - The naming system is the following: `dataset.plmod.signaltype_variabletype.weights.rds`:   
    `signaltype`: bin (binary) / cont (continuous) ;  
    `variabletype`: bs (number of binding sites) / sc (TS scores) ;  
    `weights`: w (w/ logCPM weights) / nw (w/out logCPM weights)    

The files are stored on a local server under `/mnt/schratt/enrichMiR_benchmark/results_plmod/`.

