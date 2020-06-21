for f in $(find *.fastq.gz | sed s/_[12].fastq.gz// | sort -u); do
echo $f
/conda/pkgs/star-2.6.1b-0/bin/STAR    --runThreadN 18\
            --genomeDir /reference/Mus_musculus/GENCODE/GRCm38.p5/Annotation/Genes/genes_STARIndex\
            --readFilesIn $f\_1.fastq.gz $f\_2.fastq.gz\
            --readFilesCommand zcat\
            --outFileNamePrefix ./STAR_2/$f.\
            --quantMode GeneCounts\
            --twopassMode None\
            --chimOutType Junctions SeparateSAMold
done