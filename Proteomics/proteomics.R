library(tidyverse)
library(rstatix)
library(DT)

# Load the data
lfq <- read_tsv("Proteomics/data/FS11_combined_protein.tsv") %>%
    select(c(protein, annotation, contains("LFQ"))) %>%
    mutate(annotation = str_remove(annotation, "^MAG107_contig_[0-9]+_[0-9]+ rank: [CDE]; ")) %>%
    mutate(annotation = ifelse(grepl("rank: E", annotation), "unknown", annotation)) %>%
    filter(!grepl("sp", protein)) %>%
    pivot_longer(cols = -c(protein, annotation), names_to = "sample", values_to = "LFQ") %>%
    mutate(sample = str_remove(sample, " MaxLFQ Intensity")) %>%
    mutate(LFQ = log2(LFQ)) %>%
    mutate(LFQ = if_else(is.infinite(LFQ), NA_real_, LFQ)) %>%
    pivot_wider(names_from = sample, values_from = LFQ)

ko_ids <- read_tsv("Proteomics/data/MAG107_KO_IDs.tsv") %>%
    right_join(lfq, by = "protein")


kegg <- read_tsv("KEGG_pathways.tsv") %>%
    right_join(ko_ids, by = "KO") %>%
    select(protein, annotation, gene, pathway, KO, where(is.numeric)) %>%
    mutate(gene = str_replace(gene, "^\\w{1}", toupper)) %>%
    mutate(gene = case_when(
        grepl("Flavin-binding monooxygenase", annotation) ~ "FMO",
        protein == "MAG107_contig_30_2975" ~ "SDR",
        protein == "MAG107_contig_30_387" ~ "SadC",
        protein == "MAG107_contig_30_388" ~ "SadB",
        protein == "MAG107_contig_30_389" ~ "SadA",
        protein == "MAG107_contig_30_1117" ~ "FMO",
        protein == "MAG107_contig_30_3782" ~ "CES",
        grepl("Cytochrome P450", annotation) ~ "CYP",
        TRUE ~ gene
    )) %>%
    mutate(pathway = ifelse(!is.na(gene) & is.na(pathway), "ketone", pathway)) %>%
    relocate(annotation, .after = KO) %>%
    mutate(gene = ifelse(is.na(gene), NA_character_, paste0(gene, str_extract(protein, "_[0-9]+$"))))
