#!/bin/bash


REFDIR="/storage/Genomes/Mus_musculus/GRCm38"
REFGTF="Mus_musculus.GRCm38.84.gtf"
REFFAS="Mus_musculus.GRCm38.dna.chromosome.all.fa"
REFNAM="GRCm38"

CORES=12

module load hisat2 samtools stringtie gffcompare 


ls -1 *trimmed.gtf > stringtie_mergelist.txt

stringtie --merge -p ${CORES} \
          -G ${REFDIR}/${REFGTF} \
          -o TSC_CustomTranscriptome.gtf stringtie_mergelist.txt

cat TSC_CustomTranscriptome.gtf | grep -v "^#" | awk '$3=="transcript" {print}' | wc -l


#
# Compare the reference to the new transcriptome
#

gffcompare -r ${REFDIR}/${REFGTF} -G  \
           -o TSC_CustomTranscriptome_compared.txt TSC_CustomTranscriptome.gtf

gzip TSC_CustomTranscriptome_compared.annotated.gtf 
gzip TSC_CustomTranscriptome.gtf


echo "END_OF_SCRIPT"

