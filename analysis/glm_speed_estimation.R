library(tidyverse)
library(MASS)
library(DESeq2)
library(car)


project_path <- "/cellfile/datapublic/jkoubele/celegans_mutants/"
design_formula <- ~genotype + age

metadata <- read_tsv(paste0(project_path, "sample_annotation/sample_annotation_with_library_size.tsv"))
exon_counts <- read_tsv(paste0(project_path, "feature_counts_exons_aggregated/feature_counts_aggregated.tsv"))
intron_counts <- read_tsv(paste0(project_path, "aggregated_intronic_counts/intron_read_counts.tsv"))


metadata <- metadata |>
  mutate(genotype = fct_relevel(as.factor(genotype), "wt"),
         age = fct_relevel(as.factor(age), "young")) |>
  column_to_rownames(var = "sample_name")
metadata$sample_name <- rownames(metadata)


exon_counts <- exon_counts |>
  dplyr::select(-c("Chr", "Start", "End", "Strand", "Length")) |>
  column_to_rownames(var = "Geneid") |>
  as.matrix()


dds <- DESeqDataSetFromMatrix(countData = exon_counts,
                              colData = metadata,
                              design = design_formula)
dds <- DESeq(dds)
resultsNames(dds)

deseq_results <- list()
for (result_name in resultsNames(dds)) {
  if (result_name == 'Intercept') {
    next
  }
  deseq_results[[result_name]] <- lfcShrink(dds, coef = result_name, type = "apeglm")
}


intron_counts <- intron_counts |>
  dplyr::select(-c("intron_name", "intron_name", "intron_number", "chromosome", "start", "end", "strand", "length")) |>
  group_by(gene) |>
  summarize(across(everything(), sum)) |>
  column_to_rownames("gene") |>
  filter(!if_all(everything(), ~. == 0))

# intron_counts <- head(intron_counts, 1000)

# TODO: move to DeSeq2 shrinkage for-loop
glm_intron_results <- list()
for (result_name in resultsNames(dds)) {
  if (result_name == 'Intercept') {
    next
  }
  glm_intron_results[[result_name]] <- tibble(
    gene = character(0),
    estimate = numeric(0),
    deseq_lfc = numeric(0),
    estimate_corrected_by_expression = numeric(0),
    p_value_against_zero = numeric(0),
    p_value_against_deseq_lfc = numeric(0)
  )
}


for (i in seq_len(nrow(intron_counts))) {
  if (i %% 100 == 0) {
    print(i)
  }
  row <- intron_counts[i, , drop = FALSE]
  gene_name <- rownames(row)
  
  design_matrix <- row |>
    pivot_longer(
      cols = everything(),
      names_to = "sample_name",
      values_to = "read_count") |>
    left_join(metadata, by = 'sample_name')
  
  # model <- glm.nb(formula=update(design_formula, read_count ~ . + offset(log(library_size))),
  #     data=design_matrix)
  model <- tryCatch({
    glm.nb(formula = update(design_formula, read_count ~ . + offset(log(library_size))),
           data = design_matrix)
  }, error = function(e) {
    NULL
  })
  if (is.null(model)) {
    next
  }
  model_summary <- summary(model)
  
  for (result_name in names(deseq_results)) {
    result_name_split <- strsplit(result_name, "_")[[1]]
    glm_coefficient_name <- paste0(result_name_split[1], result_name_split[2])
    coefficient_estimate <- model_summary$coefficients[glm_coefficient_name, "Estimate"]
    
    # Convert DeSeq2 L2FC to LFC with natural log:  
    deseq_lfc <- tryCatch({ deseq_results[[result_name]]['log2FoldChange'][[gene_name, 1]] / log2(exp(1)) },
                          error = function(e) {
                            NULL
                          })
    if (is.na(deseq_lfc)) {
      next
    }
    
    hypothesis_test <- linearHypothesis(model, paste(glm_coefficient_name, "=", format(deseq_lfc, scientific = FALSE)))
    
    glm_results <- tibble(gene = gene_name,
                          estimate = coefficient_estimate,
                          deseq_lfc = deseq_lfc,
                          estimate_corrected_by_expression = coefficient_estimate - deseq_lfc,
                          p_value_against_zero = model_summary$coefficients[glm_coefficient_name, "Pr(>|z|)"],
                          p_value_against_deseq_lfc = hypothesis_test$"Pr(>Chisq)"[[2]])
    
    glm_intron_results[[result_name]] <- bind_rows(glm_intron_results[[result_name]], glm_results)
  }
  
}

for (constrast_name in names(glm_intron_results)) {
  glm_results_for_contrast <- glm_intron_results[[constrast_name]]
  glm_results_for_contrast$padjust <- p.adjust(glm_results_for_contrast$p_value_against_deseq_lfc, method = 'fdr')
  volcano_data <- as.data.frame(glm_results_for_contrast) |>
    mutate(color = case_when(
      estimate_corrected_by_expression > 0 & padjust < 0.05 ~ "Upregulated",
      estimate_corrected_by_expression < 0 & padjust < 0.05 ~ "Downregulated",
      TRUE ~ "Not Significant"
    )) |>
    mutate(padjust = pmax(padjust, 1e-20),
           estimate_corrected_by_expression = pmin(pmax(estimate_corrected_by_expression, -10), 10))
  
  print(constrast_name)
  print(dplyr::count(volcano_data, color))
  print('-----------------------------------------')
  
  plot <- ggplot(volcano_data, aes(x = estimate_corrected_by_expression, y = -log10(padjust), color = color)) +
    geom_point(alpha = 0.5) +  # Set point transparency
    scale_color_manual(values = c("Not Significant" = "grey", "Upregulated" = "red", "Downregulated" = "blue")) +
    labs(title = paste("Differential Intronic Expression of", constrast_name),
         x = "Log Fold Change",
         y = "-Log10 FDR",
         color = 'Change') +
    geom_hline(yintercept = -log10(0.05), linetype = "dashed", color = "orange") +  # Add significance line
    theme_minimal() +
    theme(plot.title = element_text(hjust = 0.5))  # Center the title
  print(plot)
}