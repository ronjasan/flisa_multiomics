#!/bin/bash

mamba activate dram

samples=("S1C30" "S1C40" "S1LMWPE" "S2C30" "S2C40" "S2LMWPE")

# Rename .fna files based on the new names generated in 2_rename.R
for sample in "${samples[@]}"; do
    names_file="Metagenomics/data/new_names/${sample}_names.tsv"
    # Rename the .fna files
    while IFS=$'\t' read -r bin_id mag; do
        # Skip the header line
        if [ "$bin_id" != "Bin Id" ]; then
            # Path to the .fna file
            fna_file="DNAseq/aviary/${sample}/bins/final_bins/${mag}.fna"
            
            # Check if the .fna file exists
            if [ -f "$fna_file" ]; then
                # Rename the .fna file
                mv "$fna_file" "DNAseq/aviary/${sample}/bins/final_bins/${bin_id}.fna"
            fi
        fi
    done < $names_file
done


# Run DRAM annotation and distillation on renamed bins
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'DRAM.py annotate \
-i "DNAseq/aviary/{}/bins/final_bins/*.fna" -o "DNAseq/DRAM/{}" \
&& DRAM.py distill -i "DNAseq/DRAM/{}/annotations.tsv" -o "DNAseq/DRAM/{}/genome_summaries" \
--trna_path "DNAseq/DRAM/{}/trnas.tsv" --rrna_path "DNAseq/DRAM/{}/rrnas.tsv"'
