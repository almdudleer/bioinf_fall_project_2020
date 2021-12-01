#!/usr/bin/env Rscript
library(tidyverse)

args <- commandArgs(trailingOnly = TRUE)

if (length(args) != 3) {
  stop("Exactly three arguments must be supplied (input dir, file pattern, output file)", call. = FALSE)
}

dfs <- list.files(args[1], pattern = args[2]) %>%
  map(function(x) file.path(args[1], x)) %>%
  map(function(x) read_tsv(x, show_col_types = FALSE))

df <- dfs %>% reduce(left_join, by = "miRNA")

write.table(df, file = args[3], quote = FALSE, sep = '\t', row.names = FALSE)

