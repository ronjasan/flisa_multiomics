# Aviary
```
mamba activate aviary

samples=("S1C30" "S1C40" "S1PELW" "S2C30" "S2C40" "S2PELW")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'aviary recover --longreads DNAseq/longreads/{}.fastq.gz --long_read_type ont_hq --medaka-model r1041_e82_400bps_hac_g632 -o DNAseq/aviary/{}'
```

# DRAM
```
mamba activate dram
samples=("S1C30" "S1C40" "S1PELW" "S2C30" "S2C40" "S2PELW")

printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'DRAM.py annotate -i "DNAseq/aviary/{}/bins/final_bins/*.fna" -o "DNAseq/DRAM/{}" && DRAM.py distill -i "DNAseq/DRAM/{}/annotations.tsv" -o "DNAseq/DRAM/{}/genome_summaries" --trna_path "DNAseq/DRAM/{}/trnas.tsv" --rrna_path "DNAseq/DRAM/{}/rrnas.tsv"'
```

