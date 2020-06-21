for f in ./raw/*.fastq.gz; do
echo $f
/conda/pkgs/star-2.6.1b-0/bin/STAR    --runThreadN 18\
            --genomeDir /reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes_STARIndex\
            --readFilesIn $f\
            --readFilesCommand zcat\
            --outFileNamePrefix ./STAR_2/`basename $f .fastq.gz`\
            --quantMode GeneCounts\
            --twopassMode None\
            --chimOutType Junctions SeparateSAMold
done