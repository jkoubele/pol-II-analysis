library(argparse)
library(tidyverse)
library(stats)
library(rtracklayer)
library(rjson)
library(progress)

parser <- ArgumentParser()
parser$add_argument("--input_folder", help = "Folder containing coverage files.",
                    default = "/cellfile/datapublic/jkoubele/fli_dr_mice/coverage/no003-1_OA3",
                    required = FALSE)
parser$add_argument("--introns_file", help = "File with introns in .bed format.",
                    default = "/cellfile/datapublic/jkoubele/reference_genomes/GRCm39/introns.bed",
                    required = FALSE)
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.", required = FALSE,
                    default = "/cellfile/datapublic/jkoubele/fli_dr_mice/intron_slopes/no003-1_OA3")
args <- parser$parse_args()

fit_model_on_intron <- function(intron_row, strand_coverages, padding = 20, min_length = 50) {
  start <- intron_row$start
  end <- intron_row$end
  strand <- intron_row$strand
  chromosome <- intron_row$chromosome
  length <- intron_row$length
  if (intron_length <= min_length || intron_length <= 2 * padding) {
    return(tibble(intercept = NA, slope = NA, avg_coverage = NA, r_squared = NA))
  }
  intron_coverage_rle <- stats::window(strand_coverages[[strand]][[chromosome]],
                                       start = start + padding,
                                       end = end - padding)
  intron_coverage <- as.numeric(intron_coverage_rle)
  avg_coverage <- mean(intron_coverage)

  if (avg_coverage == 0) {
    return(tibble(intercept = NA, slope = NA, avg_coverage = avg_coverage, r_squared = NA))
  }

  x <- if (strand == '+') 1:length(intron_coverage) else length(intron_coverage):1

  model <- lm(intron_coverage ~ x)
  model_summary <- summary(model)

  return(tibble(intercept = model_summary$coefficients["(Intercept)", 'Estimate'],
                slope = model_summary$coefficients["x", 'Estimate'],
                avg_coverage = avg_coverage,
                r_squared = model_summary$r.squared))
}

input_folder <- args$input_folder
output_folder <- args$output_folder
read_counts <- fromJSON(file = paste0(input_folder, '/read_counts.json'))
library_size <- read_counts$selected_reads_forward + read_counts$selected_reads_reverse

introns <- read.table(args$introns_file,
                      col.names = c("chromosome", "start", "end", "name", "score", "strand"))
introns$chromosome <- as.character(introns$chromosome)
introns$strand <- as.character(introns$strand)
introns$length <- introns$end - introns$start

files_specification_read_pairs <- list(coverage_file_forward = paste0(input_folder, '/coverage_forward_pairs.bedGraph'),
                                       coverage_file_reverse = paste0(input_folder, '/coverage_reverse_pairs.bedGraph'),
                                       output_file_name = paste0(output_folder, '/slopes_read_pairs.tsv'))

files_specification_nascent_introns <- list(coverage_file_forward = paste0(input_folder, '/coverage_forward_nascent_introns.bedGraph'),
                                            coverage_file_reverse = paste0(input_folder, '/coverage_reverse_nascent_introns.bedGraph'),
                                            output_file_name = paste0(output_folder, '/slopes_nascent_introns.tsv'))

for (files_specification in list(files_specification_read_pairs, files_specification_nascent_introns)) {
  bed_graph_forward <- import.bedGraph(files_specification$coverage_file_forward)
  bed_graph_reverse <- import.bedGraph(files_specification$coverage_file_reverse)

  strand_coverages_of_read_pairs <- list('+' = coverage(bed_graph_forward, weight = bed_graph_forward$score / library_size * 1e6),
                                         '-' = coverage(bed_graph_reverse, weight = bed_graph_reverse$score / library_size * 1e6))

  compute_with_coverage <- partial(fit_model_on_intron, strand_coverages = strand_coverages_of_read_pairs)
  pb <- progress_bar$new(
    format = "  Computing intron slopes: [:bar] :percent in :elapsed",
    total = nrow(introns),
    clear = TRUE
  )

  compute_with_progress_bar <- function(...) {
    pb$tick()
    compute_with_coverage(...)
  }

  slopes_df <- bind_rows(pmap(introns, ~compute_with_progress_bar(list(...))))
  concat_df <- cbind(introns[c("chromosome", "start", "end", "strand", "length")], slopes_df)
  write.table(concat_df, files_specification$output_file_name, sep = "\t", row.names = F)
}
