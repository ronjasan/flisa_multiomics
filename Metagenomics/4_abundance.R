library(tidyverse)
library(MoMAColors)
library(DT)
library(data.table)
library(cowplot)
library(readxl)

samples <- c("S1C30", "S2C30", "S1C40", "S2C40", "S1PELW", "S2PELW")

for (sample in samples) {
    bins <- read_delim(paste0("/glittertind/home/ronjasan/spaceface/flisa/DNAseq/aviary/", sample, "/bins/bin_info.tsv")) %>%
        select(c(`Bin Id`, paste0(sample, ".fastq.gz Relative Abundance (%)"), classification, `Completeness (CheckM2)`, `Contamination (CheckM2)`)) %>%
        separate(classification, c("domain", "phylum", "class", "order", "family", "genus", "species"), sep = ";") %>%
        rename(MAG = `Bin Id`, abundance = paste0(sample, ".fastq.gz Relative Abundance (%)"), completeness = `Completeness (CheckM2)`, contamination = `Contamination (CheckM2)`) %>%
        mutate(across(c(3:9), ~ gsub(pattern = "^([a-z])(_{2})", replacement = "", .))) %>%
        mutate(across(where(is.character), ~ na_if(., ""))) %>%
        mutate(across(c(3:9), ~ gsub(pattern = "[\\s_][A-Z]", replacement = "", .))) %>%
        mutate(abundance = as.numeric(abundance)) %>%
        mutate(species = ifelse(grepl("[A-Z]$", species), sub(".$", "", species), species)) %>%
        mutate(community = paste0(sample)) %>%
        mutate(MIMAG = case_when(completeness > 90 & contamination < 5 ~ "High quality", completeness >= 50 & contamination < 10 ~ "Medium quality", TRUE ~ "Low quality"))
    write_tsv(bins, paste0("Metagenomics/data/aviary/", sample, "_bins.tsv"))
    assign(paste0(sample, "_bins"), bins, envir = .GlobalEnv)
}

all_bins <- rbind(S1PELW_bins, S2PELW_bins, S1C30_bins, S2C30_bins, S1C40_bins, S2C40_bins) %>%
    mutate(community = fct_relevel(community, c("S1PELW", "S2PELW", "S1C30", "S2C30", "S1C40", "S2C40"))) %>%
    relocate(abundance, .after = species) %>%
    filter(MIMAG != "Low quality")

write_tsv(all_bins, "Metagenomics/data/aviary/aviary_summary.tsv")
