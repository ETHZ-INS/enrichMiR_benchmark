# enrichMiR Benchmark

This repository contains scripts and files involved in benchmarking the enrichMiR tool. Dataset generation is contained in either `Rmd` or `sh` files. Datasets are stored as `SummarizedExperiment` (SE) objects inside `rds` files.

- [scripts for read mapping & quantification:](processing/)  
`processing/ * .STAR.sh`  
`processing/ * .datagen.raw.SE.sh`

- scripts for DEA:  
`* .datagen.DEA.SE.Rmd`  

- [SE datasets:](data/)  
`data/ * .DEA.SE.rds`  
`data/raw/ * .raw.SE.rds`  

- [TargetScan annotations & tissue-specific miRNA expression profiles:](data/)  
`data/`  

- [enrichMiR results & benchmarking scores:](results/)  
`results/ * .enrichMiR.rds`  
`results/ * .benchmark.rds`  

- a combined dataframe of all benchmark scores can be found under:  
[results/combined.benchmark.rds](results/combined.benchmark.rds)  

- a project report can be found under:  
[tgermade_thesis/tgermade_thesis.pdf](tgermade_thesis/tgermade_thesis.pdf)
