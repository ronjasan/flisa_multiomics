library(tidyverse)
library(rstatix)
samples <- c("S1C30", "S2C30", "S1C40", "S2C40", "S1PELW", "S2PELW")

# Read and concatenate the .tsv files
for (sample in samples) {
    directories <- list.dirs(path = paste0("/glittertind/home/ronjasan/spaceface/flisa/RNAseq/kallisto/kallisto_results/", sample, "/."), recursive = FALSE)
    all_kallisto <- purrr::map_df(directories, function(dir) {
        file_path <- file.path(dir, "abundance.tsv")
        if (file.exists(file_path)) {
            data <- read_tsv(file_path)
            data$directory <- basename(dir)
            return(data)
        }
    })
    kallisto <- all_kallisto %>%
        select(target_id, tpm, directory) %>%
        pivot_wider(names_from = directory, values_from = tpm) %>%
        rename(transcript_id = target_id) %>%
        select(-contains("Glc"), contains("Glc"))
    assign(paste0(sample, "_kallisto"), kallisto)
    write_tsv(kallisto, paste0("Metatranscriptomics/data/kallisto/", sample, "_kallisto.tsv"))
}


# Log-transform, filter and impute missing values in the data
source("impute_normal.R")

for (sample in samples) {
    imputed <- read_tsv(paste0("Metatranscriptomics/data/kallisto/", sample, "_kallisto.tsv")) %>%
        filter(
            rowSums(across(3:5, ~ . == 0)) <= 1 & # First three numeric columns
                rowSums(across((ncol(.) - 2):ncol(.), ~ . == 0)) <= 1 # Last three numeric columns
        ) %>%
        mutate(across(where(is.numeric), ~ log2(.))) %>%
        mutate(across(where(is.numeric), ~ impute_normal((.)))) %>%
        mutate_at(vars(c(2:ncol(.))), ~ as.numeric(.))
    assign(paste0(sample, "_imputed"), imputed)
}

# Perform t-tests on the imputed data
for (sample in samples) {
    imputed <- get(paste0(sample, "_imputed"))
    ttest <- imputed %>%
        pivot_longer(cols = -transcript_id, names_to = "Sample", values_to = "tpm") %>%
        mutate(Site = str_extract(Sample, "S[0-9]+")) %>%
        mutate(Condition = case_when(
            str_detect(Sample, "C30") ~ "C30",
            str_detect(Sample, "C40") ~ "C40",
            str_detect(Sample, "PE") ~ "PELW",
            str_detect(Sample, "Glc") ~ "Glc",
        )) %>%
        mutate(replicate = str_extract(Sample, "[0-9]$")) %>%
        group_by(transcript_id) %>%
        t_test( ## Use t_test to get accurate p-values
            tpm ~ Condition,
            p.adjust.method = "fdr",
            paired = TRUE ## Adjust as needed based on samples to be compared
        )

    write_tsv(ttest, paste0("Metatranscriptomics/data/ttest/", sample, "_ttest.tsv"))
}
