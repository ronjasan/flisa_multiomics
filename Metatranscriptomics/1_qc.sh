#!/bin/bash

# Quality check with FastQC and trim with fastp
samples=("S1C30" "S1C40" "S1LMWPE" "S1Glc" "S2C30" "S2C40" "S2LMWPE" "S2Glc")
reps=("1" "2" "3")

# FastQC
for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      fastqc RNAseq/reads/${sample}_${rep}_R1.fq.gz -o /RNAseq/fastqc
   done
done


# fastp
samples=("S1C30" "S1C40" "S1LMWPE" "S1Glc" "S2C30" "S2C40" "S2LMWPE" "S2Glc")
reps=("1" "2" "3")

for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      input1="RNAseq/reads/${sample}_${rep}_R1.fq.gz"
      input2="RNAseq/reads/${sample}_${rep}_R2.fq.gz"
      fastp -q -i ${input1} -I ${input2} \
      -o "RNAseq/fastp/${sample}/${sample}_${rep}_R1.fq.gz" \
      -O "RNAseq/fastp/${sample}/${sample}_${rep}_R2.fq.gz"
   done
done
