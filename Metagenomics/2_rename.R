library(tidyverse)

samples <- c("S1C30", "S2C30", "S1C40", "S2C40", "S1PELW", "S2PELW")

# Read bin info
for (sample in samples) {
    bins <- read_tsv(paste0("DNAseq/aviary/", sample, "/bins/bin_info.tsv")) %>%
        select(`Bin Id`) %>%
        mutate(sample = paste0(sample))
    assign(paste0(sample, "_bins"), bins, envir = .GlobalEnv)
}

all_bins <- rbind(S1C30_bins, S2C30_bins, S1C40_bins, S2C40_bins, S1PELW_bins, S2PELW_bins)

# Create new MAG names
new_names <- all_bins %>%
    select(c(`Bin Id`, sample)) %>%
    group_by(sample) %>%
    mutate(MAG = paste0("MAG", str_replace_all(sample, "\\D", ""), sprintf("%02d", row_number()))) %>%
    ungroup() %>%
    mutate(MAGs = str_replace_all(MAG, " ", "_")) %>%
    select(c(`Bin Id`, MAG, sample))

new_names_split <- split(new_names, new_names$sample)

# Write each subset to a separate TSV file
lapply(names(new_names_split), function(x) {
    write_tsv(new_names_split[[x]][, c("Bin Id", "MAG")], paste0("Metagenomics/data/new_names/", x, "_names.tsv"))
})
