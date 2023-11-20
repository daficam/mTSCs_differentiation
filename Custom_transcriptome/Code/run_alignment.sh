#!/bin/bash


REFDIR="/storage/Genomes/Mus_musculus/GRCm38"
REFGTF="Mus_musculus.GRCm38.84.gtf"
REFFAS="Mus_musculus.GRCm38.dna.chromosome.all.fa"
REFNAM="GRCm38"

CORES=12

module load hisat2 samtools stringtie


for i in *_trimmed.fq.gz;
do

  echo "hisat2:" ${i}
  hisat2 -p ${CORES} --dta -x ${REFDIR}/${REFNAM} -U ${i} -S ${i/.fq.gz/.sam} --summary-file ${i/.fq.gz/.hisat2.summary.txt}

  echo "samtools:" ${i}
  samtools sort -@ ${CORES} -o ${i/.fq.gz/.bam} ${i/.fq.gz/.sam}
  rm ${i/.fq.gz/.sam}

  echo "stringtie:" ${i}
  stringtie ${i/.fq.gz/.bam} -l ${i/_trimmed.fq.gz/} -p ${CORES} -G ${REFDIR}/${REFGTF} -o ${i/.fq.gz/.gtf}

done

echo "END_OF_SCRIPT"

