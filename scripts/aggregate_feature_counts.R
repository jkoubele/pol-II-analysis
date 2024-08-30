library(tidyverse)
library(argparse)


parser <- ArgumentParser()
parser$add_argument("--input_folder", help = "Folder containing output of featureCounts", required = TRUE)
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.", required = TRUE)
args <- parser$parse_args()

sample_names <- list.dirs(path = args$input_folder, full.names = FALSE, recursive = FALSE)
template_tibble <- read_tsv(paste0(args$input_folder, '/', sample_names[1], '/feature_counts.tsv'), skip = 1)[1:6]

load_sample_counts <- function(sample_name) {
  return(read_tsv(paste0(args$input_folder, '/', sample_name, '/feature_counts.tsv'), skip = 1)[[7]])
}

sample_counts <- sample_names |> map(load_sample_counts)
sample_counts <- set_names(sample_counts, sample_names)

output_tibble <- bind_cols(template_tibble, tibble(!!!sample_counts))
write_tsv(output_tibble, paste0(args$output_folder, '/feature_counts_aggregated.tsv'))
