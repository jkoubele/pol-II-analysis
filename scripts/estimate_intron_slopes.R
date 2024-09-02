library(argparse)
library(tidyverse)
library(stats)
library(rtracklayer)
library(rjson)
library(progress)

parser <- ArgumentParser()
parser$add_argument("--input_folder", help = "Folder containing coverage files.", required = TRUE, 
                    default='/cellfile/datapublic/jkoubele/fli_dr_mice/coverage/no004-0_OD1')
parser$add_argument("--introns_file", help = "File with introns in .bed format.", required = TRUE,
                    default='/cellfile/datapublic/jkoubele/reference_genomes/GRCm39/introns.bed')
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.", required = TRUE,
                    default='/cellfile/datapublic/jkoubele/fli_dr_mice/intron_slopes/no004-0_OD1')
args <- parser$parse_args()


compute_slope_from_definition <- function(intron_row,
                                          strand_coverages,
                                          library_size,
                                          padding = 20,
                                          min_length = 50){
  start <- intron_row$start
  end <- intron_row$end
  strand <- intron_row$strand
  chromosome <- intron_row$chromosome
  length <- intron_row$length
  if (length <= min_length || length <= 2 * padding) {
    return(tibble(slope=NA,
                  num_polymerases_per_million_reads=NA,
                  coverage_5_prime=NA,
                  coverage_5_prime_without_padding=NA,
                  coverage_3_prime=NA,
                  padding=NA,
                  length_without_padding=NA))
  }
  if(strand == '+'){
    coverage_5_prime <- as.numeric(strand_coverages[[strand]][[chromosome]][start+padding])
    coverage_5_prime_without_padding <- as.numeric(strand_coverages[[strand]][[chromosome]][start+1])
    coverage_3_prime <- as.numeric(strand_coverages[[strand]][[chromosome]][end-padding])
  }
  else if (strand=='-'){
    coverage_5_prime <- as.numeric(strand_coverages[[strand]][[chromosome]][end-padding])
    coverage_5_prime_without_padding <- as.numeric(strand_coverages[[strand]][[chromosome]][end-1])
    coverage_3_prime <- as.numeric(strand_coverages[[strand]][[chromosome]][start+padding])
  }
  
  return(tibble(slope = -(coverage_5_prime-coverage_3_prime) / (length - 2*padding) / library_size * 1e6,
                num_polymerases_per_million_reads=(coverage_5_prime-coverage_3_prime) / library_size * 1e6,
                coverage_5_prime=coverage_5_prime,
                coverage_5_prime_without_padding=coverage_5_prime_without_padding,
                coverage_3_prime=coverage_3_prime,
                padding=padding,
                length_without_padding=(length - 2*padding)))
  
}

fit_model_on_intron <- function(intron_row, strand_coverages, padding = 20, min_length = 50, apply_cummax=FALSE) {
  start <- intron_row$start
  end <- intron_row$end
  strand <- intron_row$strand
  chromosome <- intron_row$chromosome
  length <- intron_row$length
  if (length <= min_length || length <= 2 * padding) {
    return(tibble(intercept = NA, slope = NA, avg_coverage = NA, r_squared = NA))
  }
  intron_coverage_rle <- stats::window(strand_coverages[[strand]][[chromosome]],
                                       start = start + padding,
                                       end = end - padding)
  intron_coverage <- as.numeric(intron_coverage_rle)
  if (apply_cummax){
    if (strand=='+'){
      intron_coverage <- rev(cummax(rev(intron_coverage)))
    }
    else if (strand=='-'){
      intron_coverage <- cummax(intron_coverage)
    }
  }
    
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


# Compute slopes from definition:

bed_graph_forward <- import.bedGraph(paste0(input_folder, '/coverage_forward_nascent_introns.bedGraph'))
bed_graph_reverse <- import.bedGraph(paste0(input_folder, '/coverage_reverse_nascent_introns.bedGraph'))

strand_coverages_of_read_pairs <- list('+' = coverage(bed_graph_forward, weight = bed_graph_forward$score),
                                       '-' = coverage(bed_graph_reverse, weight = bed_graph_reverse$score))

compute_with_coverage <- partial(compute_slope_from_definition,
                                 strand_coverages=strand_coverages_of_read_pairs,
                                 library_size=library_size)
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
concat_df <- cbind(introns[c("chromosome", "start", "end", "strand", "name", "length")], slopes_df)
write.table(concat_df, paste0(output_folder, '/slopes_by_definition.tsv'), sep = "\t", row.names = F)

# Compute slopes by OLS:

specification_read_pairs <- list(coverage_file_forward = paste0(input_folder, '/coverage_forward_pairs.bedGraph'),
                                       coverage_file_reverse = paste0(input_folder, '/coverage_reverse_pairs.bedGraph'),
                                       output_file_name = paste0(output_folder, '/slopes_read_pairs.tsv'),
                                       apply_cummax=FALSE)

specification_read_pairs_cummax <- list(coverage_file_forward = paste0(input_folder, '/coverage_forward_pairs.bedGraph'),
                                       coverage_file_reverse = paste0(input_folder, '/coverage_reverse_pairs.bedGraph'),
                                       output_file_name = paste0(output_folder, '/slopes_read_pairs_cummax.tsv'),
                                       apply_cummax=TRUE)

specification_nascent_introns <- list(coverage_file_forward = paste0(input_folder, '/coverage_forward_nascent_introns.bedGraph'),
                                            coverage_file_reverse = paste0(input_folder, '/coverage_reverse_nascent_introns.bedGraph'),
                                            output_file_name = paste0(output_folder, '/slopes_nascent_introns.tsv'),
                                            apply_cummax=FALSE)

for (specification in list(specification_read_pairs, specification_read_pairs_cummax, specification_nascent_introns)) {
  bed_graph_forward <- import.bedGraph(specification$coverage_file_forward)
  bed_graph_reverse <- import.bedGraph(specification$coverage_file_reverse)

  strand_coverages_of_read_pairs <- list('+' = coverage(bed_graph_forward, weight = bed_graph_forward$score / library_size * 1e6),
                                         '-' = coverage(bed_graph_reverse, weight = bed_graph_reverse$score / library_size * 1e6))

  compute_with_coverage <- partial(fit_model_on_intron,
                                   strand_coverages = strand_coverages_of_read_pairs,
                                   apply_cummax=specification$apply_cummax)
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
  concat_df <- cbind(introns[c("chromosome", "start", "end", "strand", "name", "length")], slopes_df)
  write.table(concat_df, specification$output_file_name, sep = "\t", row.names = F)
}
