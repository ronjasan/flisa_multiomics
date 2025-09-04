# Aviary
```
mamba activate aviary

samples=("S1C30" "S1C40" "S1PELW" "S2C30" "S2C40" "S2PELW")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'aviary recover --longreads longreads/{}.fastq.gz --long_read_type ont_hq --medaka-model r1041_e82_400bps_hac_g632 -o aviary/{}'
```

# DRAM
```
mamba activate dram
samples=("S1C30" "S1C40" "S1PELW" "S2C30" "S2C40" "S2PELW")

printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'DRAM.py annotate -i "aviary/{}/bins/final_bins/*.fna" -o "DRAM/{}" && DRAM.py distill -i "DRAM/{}/annotations.tsv" -o "DRAM/{}/genome_summaries" --trna_path "DRAM/{}/trnas.tsv" --rrna_path "DRAM/{}/rrnas.tsv"'
```

