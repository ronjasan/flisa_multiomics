library(tidyverse)
library(rstatix)
library(DT)

# Load the data
lfq <- read_tsv("Proteomics/data/FS11_combined_protein.tsv") %>%
    select(c(Protein, Annotation, contains("LFQ"))) %>%
    mutate(Annotation = gsub("^MAG107_contig_[0-9]+_[0-9]+ rank: [CDE]; ", "", Annotation)) %>%
    mutate(Annotation = ifelse(grepl("rank: E", Annotation), "Unknown", Annotation)) %>%
    filter(!grepl("sp", Protein)) %>%
    pivot_longer(cols = -c(Protein, Annotation), names_to = "Sample", values_to = "LFQ") %>%
    mutate(Sample = str_remove(Sample, " MaxLFQ Intensity")) %>%
    mutate(LFQ = log2(LFQ)) %>%
    mutate(LFQ = if_else(is.infinite(LFQ), NA_real_, LFQ)) %>%
    pivot_wider(names_from = Sample, values_from = LFQ)

KO_IDs <- read_tsv("Proteomics/data/MAG107_KO_IDs.tsv") %>%
    right_join(lfq, by = "Protein")


KEGG <- read_tsv("Proteomics/data/KEGG_pathways.tsv") %>%
    right_join(KO_IDs, by = "KO") %>%
    select(Protein, Annotation, Gene, Pathway, KO, where(is.numeric)) %>%
    mutate(Gene = str_replace(Gene, "^\\w{1}", toupper)) %>%
    mutate(Gene = case_when(
        grepl("Flavin-binding monooxygenase", Annotation) ~ "FMO",
        Protein == "MAG107_contig_30_2975" ~ "SDR",
        Protein == "MAG107_contig_30_387" ~ "SadC",
        Protein == "MAG107_contig_30_388" ~ "SadB",
        Protein == "MAG107_contig_30_389" ~ "SadA",
        Protein == "MAG107_contig_30_1117" ~ "FMO",
        Protein == "MAG107_contig_30_3782" ~ "CES",
        grepl("Cytochrome P450", Annotation) ~ "CYP",
        TRUE ~ Gene
    )) %>%
    mutate(Pathway = ifelse(!is.na(Gene) & is.na(Pathway), "ketone", Pathway)) %>%
    relocate(Annotation, .after = KO) %>%
    arrange(Pathway) %>%
    mutate(across(where(is.numeric), as.character)) %>%
    mutate(across(where(is.character), ~ str_replace_all(., "\\.", ",")))

write_tsv(KEGG, "/glittertind/home/ronjasan/spaceface/flisa/Acinetobacter/proteomics/result_files/FS11_prot.tsv")
## Heatmap

heat <- KEGG %>%
    pivot_longer(cols = -c(Protein, Gene, Pathway, Annotation), names_to = "Sample", values_to = "LFQ")
