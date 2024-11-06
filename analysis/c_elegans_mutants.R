library(tidyverse)
library(DESeq2)
library(IDPmisc)


metadata <- read_tsv('/cellfile/datapublic/jkoubele/celegans_mutants/sample_annotation/sample_annotation_with_library_size.tsv')
exon_counts <- read_tsv('/cellfile/datapublic/jkoubele/celegans_mutants/feature_counts_exons_aggregated/feature_counts_aggregated.tsv')
intron_counts <- read_tsv('/cellfile/datapublic/jkoubele/celegans_mutants/aggregated_intronic_counts/intron_read_counts.tsv')

metadata <- metadata |>
  left_join(as_tibble(
      colSums(intron_counts|>dplyr::select(metadata$sample_name)), rownames = "sample_name") |>
                dplyr::rename(intron_reads = value),
    by = "sample_name") |>
  mutate(intron_fraction = intron_reads / library_size)

ggplot(metadata, aes(x = group_name, y = intron_fraction, color = group_name)) +
  geom_point(size=3) +
  labs(x = "Group Name",
       y = "Fraction of Intronic Reads",
       color='Group Name',
       title = "Intron Fraction by Group") +
  scale_y_continuous(limits = c(0, NA)) +
  theme_minimal() +
  theme(
    plot.title = element_text(hjust = 0.5, size = 16),  # Center and increase title size
    axis.title.x = element_text(size = 12),  # Increase x-axis title size
    axis.title.y = element_text(size = 12)   # Increase y-axis title size
  )

ggsave("/cellfile/datapublic/jkoubele/pol-II-analysis/analysis/intron_fraction_plot.png",width = 8, height = 6, dpi = 300, bg = "white")

intron_counts <- intron_counts |>
  mutate(across(all_of(metadata$sample_name), 
                ~ . / metadata$library_size[match(cur_column(), metadata$sample_name)] * 1e6))

group_averages <- metadata |>
  group_by(group_name) |>
  summarize(sample_names = list(sample_name), .groups = 'drop')

for (i in seq_len(nrow(group_averages))) {
  group <- group_averages$group_name[i]
  sample_names <- unlist(group_averages$sample_names[i])
  intron_counts[[group]] <- rowMeans(intron_counts[, sample_names])
}

intron_counts <- intron_counts |>
  dplyr::select(-all_of(metadata$sample_name))

intron_counts <- intron_counts |> mutate(L2FC_ama1_wt = log(ama1_old/wt_old))

ggplot(data = intron_counts, aes(x = L2FC_ama1_wt)) +
  geom_density(na.rm = TRUE, fill = "blue", alpha = 0.5) + # Set fill color and transparency
  labs(title = "Density of L2FC ama1 vs. WT", x = "L2FC ama1 vs. WT", y = "Density") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16))

ggsave("/cellfile/datapublic/jkoubele/pol-II-analysis/analysis/intron_reads_lfc_density.png",width = 8, height = 6, dpi = 300, bg = "white")

metadata_deseq <- metadata |> 
  filter(group_name %in% c('ama1_old', 'wt_old')) |>
  column_to_rownames(var = 'sample_name') |>
  mutate(group_name = fct_relevel(as.factor(group_name), "wt_old"))

exon_counts <- exon_counts |> 
  column_to_rownames(var = "Geneid") |>
  dplyr::select(rownames(metadata_deseq)) |> 
  as.matrix()

dds <- DESeqDataSetFromMatrix(countData = exon_counts,
                              colData = metadata_deseq,
                              design = ~ group_name)
dds <- DESeq(dds)
resultsNames(dds)

res <- na.omit(results(dds, name="group_name_ama1_old_vs_wt_old"))
res_shrinked <- na.omit(lfcShrink(dds, coef="group_name_ama1_old_vs_wt_old", type="apeglm"))

summary(res_shrinked)

volcano_data <- as.data.frame(res_shrinked) |>
  mutate(color = case_when(
    log2FoldChange > 0 & padj < 0.05 ~ "Upregulated",
    log2FoldChange < 0 & padj < 0.05 ~ "Downregulated",
    TRUE ~ "Not Significant"
  )) |>
  rownames_to_column(var = "gene")

ggplot(volcano_data, aes(x = log2FoldChange, y = -log10(padj), color = color)) +
  geom_point(alpha = 0.5) +  # Set point transparency
  scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
  labs(title = "DE of ama1 vs WT", 
       x = "Log2 Fold Change", 
       y = "-Log10 Adjusted P-Value",
       color='Change') +
  geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +  # Add significance line
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

ggsave("/cellfile/datapublic/jkoubele/pol-II-analysis/analysis/volcano_plot.png",width = 8, height = 6, dpi = 300, bg = "white")

intron_counts <- intron_counts |>
  left_join(volcano_data |> dplyr::select(gene, log2FoldChange), by = "gene") |>
  dplyr::rename(L2FC_ama1_wt_DE = log2FoldChange) |>
  NaRV.omit()


ggplot(intron_counts, aes(x = L2FC_ama1_wt_DE, y = L2FC_ama1_wt)) +
  geom_point(alpha = 0.45, color = 'blue') +  # Adjust transparency of points
  labs(title = "Scatter Plot of L2FC Values",
       x = "L2FC ama1 vs wt (DE)",
       y = "L2FC ama1 vs wt") +
  theme_minimal() +  # Use a minimal theme
  theme(plot.title = element_text(hjust = 0.5))  # Center the title

model <- lm(L2FC_ama1_wt ~ L2FC_ama1_wt_DE, data = intron_counts)
summary(model)

ggplot(intron_counts, aes(x = L2FC_ama1_wt_DE, y = L2FC_ama1_wt)) +
  geom_point(alpha = 0.5) + 
  geom_smooth(method = "lm", color = "blue", se = TRUE) + 
  labs(title = "L2FC of intronic reads vs. L2FC of corresponding genes",
       x = "L2FC of gene expression",
       y = "L2FC of intronic reads") +
  theme_minimal() +
  theme(plot.title = element_text(hjust = 0.5))  

ggsave("/cellfile/datapublic/jkoubele/pol-II-analysis/analysis/intron_L2FC_vs_gene_L2FC.png",width = 8, height = 6, dpi = 300, bg = "white")

intron_counts <- intron_counts |> 
  mutate(L2FC_ama1_wt_corrected = L2FC_ama1_wt - L2FC_ama1_wt_DE)

ggplot(data = intron_counts, aes(x = L2FC_ama1_wt_corrected)) +
  geom_density(na.rm = TRUE, fill = "blue", alpha = 0.5) + # Set fill color and transparency
  labs(title = "Density of corrected L2FC ama1 vs. WT", x = "Corrected L2FC ama1 vs. WT", y = "Density") +
  theme_minimal() + 
  theme(plot.title = element_text(hjust = 0.5, size = 16))

ggsave("/cellfile/datapublic/jkoubele/pol-II-analysis/analysis/intron_reads_lfc_corrected_density.png",width = 8, height = 6, dpi = 300, bg = "white")



