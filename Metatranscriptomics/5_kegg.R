library(tidyverse)

samples <- c("S1C30", "S2C30", "S1C40", "S2C40", "S1PELW", "S2PELW")

# Combine KEGG, DRAM, t-test and kallisto data for each sample
for (sample in samples) {
    mags <- read_tsv(paste0("Metagenomics/data/aviary/", sample, "_bins.tsv")) %>%
        filter(MIMAG != "Low quality")

    kallisto <- read_tsv(paste0("Metatranscriptomics/data/imputed/", sample, "_imputed.tsv")) %>%
        mutate(MAG = str_extract(transcript_id, "MAG[0-9]+")) %>%
        select(c(MAG, transcript_id, where(is.numeric))) %>%
        right_join(mags, by = "MAG") %>%
        mutate(
            mean = rowMeans(across(3:5)),
            mean_control = rowMeans(across(6:8, ))
        ) %>%
        mutate(log2FC = mean - mean_control) %>%
        select(MAG, transcript_id, log2FC)

    ttest <- read_tsv(paste0("Metatranscriptomics/data/ttest/", sample, "_ttest.tsv")) %>%
        select(transcript_id, p) %>%
        mutate(log10p = -log10(p)) %>%
        right_join(kallisto, by = "transcript_id")

    dram <- read_tsv(paste0("DNAseq/DRAM/", sample, "/annotations.tsv")) %>%
        mutate(pfam_hits = str_remove(pfam_hits, "\\[.*")) %>%
        mutate(dram_annotation = case_when(rank == "E" ~ "unknown", rank == "D" ~ pfam_hits, rank == "C" ~ kegg_hit)) %>%
        rename(transcript_id = 1, KO = ko_id) %>%
        select(transcript_id, dram_annotation, gene_position, kegg_hit, pfam_hits, KO) %>%
        right_join(ttest, by = "transcript_id")

    kegg <- read_tsv("KEGG_pathways.tsv") %>%
        select(KO, pathway, gene) %>%
        right_join(dram, by = "KO", relationship = "many-to-many") %>%
        distinct() %>%
        mutate(gene = case_when(
            is.na(kegg_hit) & grepl("Aldehyde dehydrogenase family", dram_annotation) ~ "ALDH",
            is.na(kegg_hit) & grepl("Rubredoxin", dram_annotation) ~ "rub",
            is.na(kegg_hit) & grepl("Flavin-binding monooxygenase-like", dram_annotation) ~ "FMO",
            is.na(kegg_hit) & grepl("Alcohol dehydrogenase GroES-like domain", dram_annotation) ~ "ADH",
            is.na(kegg_hit) & grepl("Cytochrome P450", dram_annotation) ~ "CYP",
            TRUE ~ gene
        )) %>%
        mutate(pathway = case_when(
            is.na(pathway) & grepl("ALDH", gene) ~ "alkane",
            is.na(pathway) & grepl("rub", gene) ~ "alkane",
            is.na(pathway) & grepl("FMO", gene) ~ "ketone",
            is.na(pathway) & grepl("ADH", gene) ~ "alkane",
            is.na(pathway) & grepl("CYP", gene) ~ "alkane",
            TRUE ~ pathway
        )) %>%
        mutate(gene_name = case_when(
            !is.na(gene) ~ paste0(MAG, "_", gene, "_", gene_position),
            TRUE ~ NA_character_
        )) %>%
        select(transcript_id, pathway, gene, gene_name, dram_annotation, contains(c(paste0(sample), "Glc")), log2FC, log10p)
    assign(paste0(sample, "_kegg"), kegg)
}
