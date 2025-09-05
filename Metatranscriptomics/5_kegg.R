library(tidyverse)
library(MoMAColors)
library(rstatix)
library(ggrepel)
library(cowplot)
library(DT)
library(ggnewscale)

samples <- c("S1C30", "S2C30", "S1C40", "S2C40", "S1PELW", "S2PELW")

for (sample in samples) {
    mags <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/DNAseq/R/MAG_info/", sample, "_summary.tsv")) %>%
        filter(MIMAG != "Low quality")

    kallisto <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/imp_data/", sample, "_imputed.tsv")) %>%
        mutate(MAG = str_extract(transcript_id, "MAG[0-9]+")) %>%
        select(c(MAG, transcript_id, where(is.numeric))) %>%
        right_join(mags, by = "MAG") %>%
        mutate(
            mean = rowMeans(across(3:5)),
            mean_control = rowMeans(across(6:8, ))
        ) %>%
        mutate(log2FC = mean - mean_control) %>%
        select(MAG, transcript_id, log2FC)

    ttest <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/t_test_data/", sample, "_ttest.tsv")) %>%
        select(transcript_id, p) %>%
        mutate(log10p = -log10(p)) %>%
        right_join(kallisto, by = "transcript_id")

    dram <- read_tsv(paste0("/spaceface/projects/enzyclic/flisa/DNAseq/DRAM/", sample, "/annotations.tsv")) %>%
        mutate(pfam_hits = str_remove(pfam_hits, "\\[.*")) %>%
        mutate(dram_annotation = case_when(rank == "E" ~ "unknown", rank == "D" ~ pfam_hits, rank == "C" ~ kegg_hit)) %>%
        rename(MAG = fasta, transcript_id = 1, KO = ko_id) %>%
        select(transcript_id, MAG, dram_annotation, gene_position, kegg_hit, pfam_hits, KO)

    relevant <- read_tsv("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/MAGs.tsv") %>%
        filter(Sample == sample)

    volcano <- read_tsv("/glittertind/home/ronjasan/spaceface/kegg_pathways/fatty_acid_degradation.tsv") %>%
        select(KO, pathway, gene) %>%
        right_join(dram, by = "KO", relationship = "many-to-many") %>%
        distinct() %>%
        right_join(relevant, by = "MAG") %>%
        mutate(gene = case_when(
            !is.na(relevance) & is.na(kegg_hit) & grepl("Acyclic terpene utilisation family protein AtuA", dram_annotation) ~ "atuA",
            !is.na(relevance) & is.na(kegg_hit) & grepl("ER-bound oxygenase mpaB/B'/Rubber oxygenase, catalytic domain", dram_annotation) ~ "rox",
            !is.na(relevance) & is.na(kegg_hit) & grepl("Aldehyde dehydrogenase family", dram_annotation) ~ "ALDH",
            !is.na(relevance) & is.na(kegg_hit) & grepl("Rubredoxin", dram_annotation) ~ "rub",
            !is.na(relevance) & is.na(kegg_hit) & grepl("Flavin-binding monooxygenase-like", dram_annotation) ~ "FMO",
            !is.na(relevance) & is.na(kegg_hit) & grepl("Alcohol dehydrogenase GroES-like domain", dram_annotation) ~ "ADH",
            !is.na(relevance) & is.na(kegg_hit) & grepl("Cytochrome P450", dram_annotation) ~ "CYP",
            TRUE ~ gene
        )) %>%
        mutate(pathway = case_when(
            is.na(pathway) & grepl("atuA", gene) ~ "other",
            is.na(pathway) & grepl("rox", gene) ~ "other",
            is.na(pathway) & grepl("ALDH", gene) ~ "alkane",
            is.na(pathway) & grepl("rub", gene) ~ "alkane",
            is.na(pathway) & grepl("FMO", gene) ~ "ketone",
            is.na(pathway) & grepl("ADH", gene) ~ "alkane",
            is.na(pathway) & grepl("CYP", gene) ~ "alkane",
            TRUE ~ pathway
        )) %>%
        mutate(gene_name = case_when(
            !is.na(gene) & !is.na(relevance) ~ paste0(MAG, "_", gene, "_", gene_position),
            TRUE ~ NA_character_
        )) %>%
        select(transcript_id, pathway, gene, gene_name) %>%
        right_join(ttest, by = "transcript_id")

    assign(paste0(sample, "_volcano"), volcano)
    # write_tsv(volcano, paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/volcano_data/", sample, "_volcano.tsv"))
}



## Volcano plots
for (sample in samples) {
    volcano <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/volcano_data/", sample, "_volcano.tsv")) %>%
        mutate(sign = case_when(
            log10p > 1.3 & log2FC > 2.5 ~ "Upregulated",
            log10p > 1.3 & log2FC < -2.5 ~ "Downregulated",
            TRUE ~ "Not signifigant"
        )) %>%
        mutate(gene_name = case_when(
            sign == "Not signifigant" ~ NA_character_,
            sign == "Downregulated" ~ NA_character_,
            TRUE ~ gene_name
        )) %>%
        mutate(pathway = case_when(
            sign == "Not signifigant" ~ NA_character_,
            sign == "Downregulated" ~ NA_character_,
            TRUE ~ pathway
        )) %>%
        mutate(pathway = as.factor(pathway))
    assign(paste0(sample, "_volcano"), volcano)
    volcano_plot <- volcano %>%
        ggplot(aes(x = log2FC, y = log10p)) +
        geom_point(size = 0.5, aes(color = sign)) +
        geom_vline(xintercept = c(-2.5, 2.5), col = "darkgray", linetype = "dashed") +
        geom_hline(yintercept = 1.3, col = "darkgray", linetype = "dashed") +
        scale_color_manual(values = c("Upregulated" = "#c41f50", "Downregulated" = "#41606d", "Not signifigant" = "grey")) +
        new_scale_color() +
        geom_label_repel(aes(label = gene_name, colour = pathway, fontface = "bold", segment.color = "black"), box.padding = 1, max.overlaps = Inf, min.segment.length = 0.1, size = 3, na.rm = TRUE) +
        scale_color_manual(values = c("ketone" = "#edc919", "alkane" = "#69a256", "other" = "#582851", "NA" = "darkgrey")) +
        labs(title = paste0(sample), x = paste0("log2(", sample, "/Glc)"), y = "-log10(p-value)") +
        theme_classic() +
        theme(
            panel.border = element_blank(),
            panel.background = element_blank(),
            axis.ticks = element_line(color = "black", linewidth = 0.75),
            axis.text = element_text(color = "black", size = 14),
            axis.ticks.length = unit(0.2, "cm"),
            axis.text.x = element_text(margin = margin(5, 0, 0, 0)),
            axis.text.y = element_text(margin = margin(0, 5, 0, 0)),
            axis.title = element_text(size = 14, face = "bold"),
            axis.title.x = element_text(margin = margin(10, 0, 0, 0)),
            axis.title.y = element_text(margin = margin(0, 10, 0, 0)),
            plot.title = element_text(size = 16, face = "bold", hjust = 0.5),
            legend.title = element_blank(),
            legend.text = element_text(size = 10, face = "bold"),
            plot.margin = margin(0.5, 0.5, 0.5, 0.5, "cm")
        )
    assign(paste0(sample, "_volcano_plot"), volcano_plot)
}


S1PELW_volcano_plot

legend_volcano <- get_legend(
    S1PELW_volcano_plot,
    return_all = TRUE
)
ggsave(legend_volcano, filename = "/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/plots/legend_volcano.pdf", width = 5, height = 2.5, units = "cm")


all_volcano <- plot_grid(
    S1PELW_volcano_plot + theme(legend.position = "none"),
    S1C30_volcano_plot + theme(legend.position = "none"),
    S1C40_volcano_plot + theme(legend.position = "none"),
    S2PELW_volcano_plot + theme(legend.position = "none"),
    S2C30_volcano_plot + theme(legend.position = "none"),
    S2C40_volcano_plot + theme(legend.position = "none"),
    ncol = 3, nrow = 2,
    labels = c("A", "B", "C", "D", "E", "F"),
    label_size = 20
)

ggsave(all_volcano, filename = "/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/plots/volcano_all.pdf", width = 45, height = 35, units = "cm")


## TEST
mags <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/DNAseq/R/MAG_info/S1PELW_summary.tsv")) %>%
    filter(MIMAG != "Low quality") %>%
    select(MAG)

kallisto <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/imp_data/S1PELW_imputed.tsv")) %>%
    mutate(MAG = str_extract(transcript_id, "MAG[0-9]+")) %>%
    relocate(MAG, .before = transcript_id) %>%
    right_join(mags, by = "MAG") %>%
    mutate(
        mean = rowMeans(across(3:5)),
        mean_control = rowMeans(across(6:8, ))
    ) %>%
    mutate(log2FC = mean - mean_control) %>%
    select(-c(mean, mean_control))

ttest <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/t_test_data/S1PELW_ttest.tsv")) %>%
    select(transcript_id, p)
dram <- read_tsv(paste0("/spaceface/projects/enzyclic/flisa/DNAseq/DRAM/S1PELW/annotations.tsv")) %>%
    mutate(pfam_hits = str_remove(pfam_hits, "\\[.*")) %>%
    mutate(dram_annotation = case_when(rank == "E" ~ "unknown", rank == "D" ~ pfam_hits, rank == "C" ~ kegg_hit)) %>%
    rename(MAG = fasta, transcript_id = 1, KO = ko_id) %>%
    select(transcript_id, dram_annotation, kegg_hit, pfam_hits, KO)

upregulated <- read_tsv("/glittertind/home/ronjasan/spaceface/kegg_pathways/fatty_acid_degradation_all.tsv") %>%
    select(KO, pathway, gene) %>%
    right_join(dram, by = "KO", relationship = "many-to-many") %>%
    distinct() %>%
    right_join(kallisto, by = "transcript_id") %>%
    right_join(ttest, by = "transcript_id") %>%
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
    select(MAG, transcript_id, pathway, gene, dram_annotation, contains(c(paste0(sample), "Glu")), log2FC, p) %>%
    filter(log2FC > 2.5 & p < 0.05) %>%
    arrange(desc(log2FC))


## All alkane + fatty acid
samples <- c("S1C30", "S2C30", "S1C40", "S2C40", "S1PELW", "S2PELW")

for (sample in samples) {
    mags <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/DNAseq/R/MAG_info/", sample, "_summary.tsv")) %>%
        filter(MIMAG != "Low quality") %>%
        select(MAG)

    kallisto <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/imp_data/", sample, "_imputed.tsv")) %>%
        mutate(MAG = str_extract(transcript_id, "MAG[0-9]+")) %>%
        relocate(MAG, .before = transcript_id) %>%
        right_join(mags, by = "MAG") %>%
        mutate(
            mean = rowMeans(across(3:5)),
            mean_control = rowMeans(across(6:8, ))
        ) %>%
        mutate(log2FC = mean - mean_control) %>%
        select(-c(mean, mean_control))

    ttest <- read_tsv(paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/t_test_data/", sample, "_ttest.tsv")) %>%
        select(transcript_id, p)
    dram <- read_tsv(paste0("/spaceface/projects/enzyclic/flisa/DNAseq/DRAM/", sample, "/annotations.tsv")) %>%
        mutate(pfam_hits = str_remove(pfam_hits, "\\[.*")) %>%
        mutate(dram_annotation = case_when(rank == "E" ~ "unknown", rank == "D" ~ pfam_hits, rank == "C" ~ kegg_hit)) %>%
        rename(MAG = fasta, transcript_id = 1, KO = ko_id) %>%
        select(transcript_id, dram_annotation, kegg_hit, pfam_hits, KO)

    upregulated <- read_tsv("/glittertind/home/ronjasan/spaceface/kegg_pathways/fatty_acid_degradation_all.tsv") %>%
        select(KO, pathway, gene) %>%
        right_join(dram, by = "KO", relationship = "many-to-many") %>%
        distinct() %>%
        right_join(kallisto, by = "transcript_id") %>%
        right_join(ttest, by = "transcript_id") %>%
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
        select(MAG, transcript_id, pathway, gene, dram_annotation, contains(c(paste0(sample), "Glc")), log2FC, p) %>%
        filter(log2FC > 2.5 & p < 0.05) %>%
        arrange(desc(log2FC)) %>%
        mutate(across(where(is.numeric), as.character)) %>%
        mutate(across(where(is.character), ~ str_replace_all(., "\\.", ",")))

    assign(paste0(sample, "_upregulated"), upregulated)
    write_tsv(upregulated, paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/R/upregulated_2/", sample, ".tsv"))
}
