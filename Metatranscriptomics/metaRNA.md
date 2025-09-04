# FastQC
```
# merge reads
samples=("S1C30" "S1C40" "S1PELW" "S1Glc" "S2C30" "S2C40" "S2PELW" "S2Glc")
reps=("1" "2" "3")

for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      zcat RNAseq/reads/${sample}/${sample}_${rep}_L{1..4}_1.fq.gz | gzip > RNAseq/reads/merged/${sample}/${sample}_${rep}_R1.fq.gz
      zcat RNAseq/reads/${sample}/${sample}_${rep}_L{1..4}_2.fq.gz | gzip > RNAseq/reads/merged/${sample}/${sample}_${rep}_R2.fq.gz
   done
done


mamba activate fastqc 

for sample in "${samples[@]}"
do
   for rep in "${reps[@]}"
   do
      fastqc RNAseq/reads/merged/${sample}/${sample}_${rep}_R1.fq.gz -o /RNAseq/fastqc
   done
done
```

# fastp
```
mamba activate fastp
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
```

# SortMeRNA
```
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
```

# Kallisto
## Index
```
samples=("S1C30" "S1C40" "S1PELW" "S2C30" "S2C40" "S2PELW")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 6 -I{} bash -c 'kallisto index -i RNAseq/kallisto/idx/{}/{}.idx DNAseq/DRAM/{}/genes.fna'
```

## Pseudo-align
```
samples=("S1C30" "S1C40" "S1PELW" "S1Glc" "S2C30" "S2C40" "S2PELW" "S2Glc")
reps=("1" "2" "3")
index=RNAseq/kallisto/idx
output=RNAseq/kallisto/results
input=RNAseq/sortmerna

printf "%s\n" "${samples[@]}" | xargs -I{} -P 8 bash -c '
  for rep in "$@"
  do
    kallisto quant -i '"$index"'/{}/{}.idx \
    -o '"$output"'/{}/{}_${rep} -b 100 \
    '"$input"'/{}_${rep}_R1.fq.gz '"$iput"'/{}${rep}_R2.fq.gz
  done
' bash "${reps[@]}"
```