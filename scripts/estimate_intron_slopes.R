library(argparse)
library(stats)
library(rtracklayer)
library(rjson)


parser <- ArgumentParser()
parser$add_argument("--input_folder", help = "Folder containing coverage files.",
                    default = "/cellfile/datapublic/jkoubele/fli_dr_mice/coverage/no003-1_OA3/",
                    required = FALSE)
parser$add_argument("--introns_file", help = "File with introns in .bed format.", 
                    default="/cellfile/datapublic/jkoubele/reference_genomes/GRCm39/introns.bed",
                    required = FALSE)
parser$add_argument("--output_folder", help = "Folder to which the result will be saved.", required = FALSE,
                    default="/cellfile/datapublic/jkoubele/fli_dr_mice/intron_slopes/no003-1_OA3")
args <- parser$parse_args()

input_folder <- args$input_folder
output_folder <- args$output_folder

read_counts <- fromJSON(file=paste0(input_folder, 'read_counts.json'))
library_size <- read_counts$selected_reads_forward + read_counts$selected_reads_reverse

introns <- read.table(args$introns_file, 
                      col.names=c("chromosome", "start", "end", "name", "score", "strand"))
introns$chromosome <- as.character(introns$chromosome)
introns$strand <- as.character(introns$strand)

coverage_file_forward <-  paste0(input_folder, 'coverage_forward_pairs.bedGraph')
coverage_file_reverse <-  paste0(input_folder, 'coverage_reverse_pairs.bedGraph')

bed_graph_forward <- import.bedGraph(coverage_file_forward)
bed_graph_reverse <- import.bedGraph(coverage_file_reverse)

strand_coverages <- list('+' = coverage(bed_graph_forward, weight = bed_graph_forward$score / library_size * 1e6), 
                         '-' = coverage(bed_graph_reverse, weight = bed_graph_reverse$score / library_size * 1e6))

slope <- c()
intercept <- c()
avg_coverage <- c()
padding <- 20
min_length <- 50

for (i in 1:nrow(introns)){
  if (i%%1000==0){
    print(i)
  }
  
  row <- introns[i,]
  intron_start <- row[['start']]
  intron_end <- row[['end']]
  intron_strand <- row[['strand']]
  intron_chromosome <- row[['chromosome']]
  if (intron_end - intron_start < min_length){
    slope <- c(slope, NA)
    intercept <- c(intercept, NA)
    avg_coverage <- c(avg_coverage, NA)
    next
  }
  
  intron_coverage_rle <- stats::window(strand_coverages[[intron_strand]][[intron_chromosome]], 
                                start=intron_start+padding,
                                end=intron_end-padding)
  intron_coverage <- as.numeric(intron_coverage_rle)
  if (mean(intron_coverage)==0){
    slope <- c(slope, NA)
    intercept <- c(intercept, NA)
    avg_coverage <- c(avg_coverage, mean(intron_coverage))
    next
  }
  
  x <-if (intron_strand=='+') 1:length(intron_coverage) else length(intron_coverage):1
  
  model <- lm(intron_coverage~x)
  model_summary <- summary(model)
  
  slope <- c(slope, model_summary$coefficients["x",'Estimate'])
  intercept <- c(intercept, model_summary$coefficients["(Intercept)",'Estimate'])
  avg_coverage <- c(avg_coverage, mean(intron_coverage))
}

result_df <- data.frame(slope = slope,
                        intercept = intercept,
                        avg_coverage = avg_coverage)

concat_df <- cbind(introns[c("chromomsome", "start", "end", "strand")], result_df)
#write.table(concat_df, paste0(output_folder, sample_name, '.tsv'), sep = "\t", row.names = F)

