# FastQC
```
mamba activate fastqc 

samples=("S1C30" "S1C40" "S1PELW" "S1Glc" "S2C30" "S2C40" "S2PELW" "S2Glc")

for sample in "${sample[@]}"
do
   fastqc RNAseq/$sample/*1.fq.gz -o /glittertind/home/ronjasan/scratch/flisa/RNAseq
done