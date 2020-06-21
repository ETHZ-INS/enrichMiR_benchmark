#!/bin/bash/

mkdir STAR

# Alignment for samples SRR10513924-SRR10513935 (paired-end RNA-seq):

for f in *_1.fastq.gz; do
/conda/pkgs/star-2.6.1b-0/bin/STAR    --runThreadN 18\
            --genomeDir /reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes_STARIndex\
            --genomeLoad LoadAndKeep\
            --readFilesIn $f ${f%_1.fastq.gz}_2.fastq.gz\
            --readFilesCommand zcat\
            --outFileNamePrefix ./STAR/${f%_1.fastq.gz}.\
            --quantMode GeneCounts\
            --twopassMode None\
	    --outFilterScoreMin 1\
            --chimOutType Junctions SeparateSAMold\
	    --chimScoreMin 1
done
