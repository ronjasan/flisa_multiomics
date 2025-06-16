# Aviary

```
mamba activate aviary
cd /glittertind/home/ronjasan/scratch/flisa/DNAseq/aviary

samples=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'aviary recover --longreads /spaceface/projects/enzyclic/flisa/DNAseq/longreads/{}.fastq.gz --long_read_type ont_hq --medaka-model r1041_e82_400bps_hac_g632 -o {}'
```

## Rename all MAGs
Files to be renamed:
{sample}/bins/final_bins/*.

Files that contain the Bin Id:
{sample}/bins/bin_info.tsv


Then rerun gtdbtk.

samples=("S1C30" "S2C30" "S1C40" "S2C40" "S1PE" "S2PE" "S3PE")
for sample in "${samples[@]}"
do
  /usr/bin/rm -r $sample/data/gtdbtk
  /usr/bin/rm -r $sample/taxonomy
done

mamba activate aviary
cd /glittertind/home/ronjasan/scratch/flisa/DNAseq/aviary

samples=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'aviary recover -w gtdbtk --assembly {}/assembly/final_contigs.fasta --longreads /spaceface/projects/enzyclic/flisa/DNAseq/longreads/{}.fastq.gz --long_read_type ont_hq --medaka-model r1041_e82_400bps_hac_g632 -o {}'

# DRAM

```
mamba activate dram
cd /glittertind/home/ronjasan/scratch/flisa/DNAseq/

samples=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")

samples=("S2C40")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 4 -I{} bash -c 'DRAM.py annotate -i "aviary/{}/bins/final_bins/*.fna" -o "DRAM/{}" && DRAM.py distill -i "DRAM/{}/annotations.tsv" -o "DRAM/{}/genome_summaries" --trna_path "DRAM/{}/trnas.tsv" --rrna_path "DRAM/{}/rrnas.tsv"'
```

# iqtree

```
dirs=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
for dir in "${dirs[@]}"
do 
  mkdir iqtree/"$dir"
done


mamba activate iqtree
cd /glittertind/home/ronjasan/scratch/flisa/DNAseq/

samples=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 7 -I{} bash -c 'iqtree -m MFP -bb 1000 -s aviary/{}/data/gtdbtk/align/gtdbtk.bac120.user_msa.fasta.gz -pre iqtree/{}/{}'

```


# HADEG

```
cd /glittertind/home/ronjasan/scratch/flisa/DNAseq/HADEG
mamba activate proteinortho

samples=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 7 -I{} bash -c 'cd {}/proteinortho && proteinortho fasta/*.faa -identity=50 -conn=0.3 -project={} --clean'
```

dirs=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
for dir in "${dirs[@]}"
do 
  /usr/bin/rm $dir.faa 
done

## SignalP
samples=("S1C30" "S1C40" "S1PE" "S2C30" "S2C40" "S2PE" "S3PE")
printf "%s\n" "${samples[@]}" | xargs -n 1 -P 7 -I{} bash -c 'signalp6 --fastafile /spaceface/projects/enzyclic/flisa/DNAseq/DRAM/{}/genes.faa --organism other --output_dir /spaceface/projects/signalp_tmp/{} --mode slow'