#!/bin/bash

# Kallisto
samples=("S1C30" "S1C40" "S1LMWPE" "S2C30" "S2C40" "S2LMWPE")
reps=("1" "2" "3")
index=RNAseq/kallisto/idx
output=RNAseq/kallisto/results
input=RNAseq/sortmerna

## Indexing
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 6 -I{} bash -c 'kallisto index \
-i RNAseq/kallisto/idx/{}/{}.idx DNAseq/DRAM/{}/genes.fna'

## Pseudoalignment and quantification
printf "%s\n" "${samples[@]}" | xargs -I{} -P 8 bash -c '
  for rep in "$@"
  do
    kallisto quant -i '"$index"'/{}/{}.idx \
    -o '"$output"'/{}/{}_${rep} -b 100 \
    '"$input"'/{}_${rep}_R1.fq.gz '"$input"'/{}_${rep}_R2.fq.gz
  done
' bash "${reps[@]}"