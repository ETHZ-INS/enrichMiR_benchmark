#!/bin/bash/

mkdir STAR

for f in *.fastq.gz; do
STAR    --runThreadN 24\
            --genomeDir /reference/Homo_sapiens/Ensembl/GRCh38.p13/Annotation/Release_98-2019-09-3/genes_STARIndex\
            --genomeLoad LoadAndKeep\
            --readFilesIn $f\
            --readFilesCommand zcat\
            --outFileNamePrefix ./STAR2/`basename $f .fastq.gz`\
            --quantMode GeneCounts\
            --twopassMode None\
            --chimOutType Junctions SeparateSAMold
done
