#!/bin/bash

# Merge raw reads, quality check with FastQC and trim with fastp
samples=("S1C30" "S1C40" "S1PELW" "S1Glc" "S2C30" "S2C40" "S2PELW" "S2Glc")
reps=("1" "2" "3")

## Merge
for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      zcat RNAseq/reads/${sample}/${sample}_${rep}_L{1..4}_1.fq.gz | gzip > RNAseq/reads/merged/${sample}/${sample}_${rep}_R1.fq.gz
      zcat RNAseq/reads/${sample}/${sample}_${rep}_L{1..4}_2.fq.gz | gzip > RNAseq/reads/merged/${sample}/${sample}_${rep}_R2.fq.gz
   done
done

# FastQC
for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      fastqc RNAseq/reads/merged/${sample}/${sample}_${rep}_R1.fq.gz -o /RNAseq/fastqc
   done
done


# fastp
samples=("S1C30" "S1C40" "S1PELW" "S1Glc" "S2C30" "S2C40" "S2PELW" "S2Glc")
reps=("1" "2" "3")

for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      input1="RNAseq/reads/merged/${sample}/${sample}_${rep}_R1.fq.gz"
      input2="RNAseq/reads/merged/${sample}/${sample}_${rep}_R2.fq.gz"
      fastp -q -i ${input1} -I ${input2} \
      -o "RNAseq/fastp/${sample}/${sample}_${rep}_R1.fq.gz" \
      -O "RNAseq/fastp/${sample}/${sample}_${rep}_R2.fq.gz"
   done
done