#!/bin/bash

# SortMeRNA
samples=("S1C30" "S1C40" "S1PELW" "S1Glc" "S2C30" "S2C40" "S2PELW" "S2Glc")
reps=("1" "2" "3")

for sample in "${samples[@]}"
do
  for rep in "${reps[@]}"
  do 
    r1="RNAseq/fastp/${sample}/${sample}_${rep}_R1.fq.gz"
    r2="RNAseq/fastp/${sample}/${sample}_${rep}_R2.fq.gz"
    if [[ -e "$r1" ]] && [[ -e "$r2" ]]
    then
      mkdir -p "RNAseq/sortmerna/${sample}_${rep}"
      touch "RNAseq/sortmerna/${sample}_${rep}/${sample}_${rep}.log"
      sortmerna --workdir "RNAseq/sortmerna/${sample}_${rep}" \
      --out2 --paired_out --fastx --other "RNAseq/sortmerna/${sample}_${rep}" --threads 24 \
      --ref databases/rRNA_databases/silva-bac-16s-id90.fasta --ref databases/rRNA_databases/silva-bac-23s-id98.fasta \
      --ref databases/rRNA_databases/silva-arc-16s-id95.fasta --ref databases/rRNA_databases/silva-arc-23s-id98.fasta \
      --reads "$r1" --reads "$r2" 
    fi
  done
done