 - First QC: ```sh batch_qc.sh -i /data/public/jkoubele/senescent_cells/FASTQ -o /data/public/jkoubele/senescent_cells/QC_before_trimming```
 - Detecting adapters: ```sh batch_detect_adapers.sh -i /data/public/jkoubele/senescent_cells/FASTQ -o /data/public/jkoubele/senescent_cells/detected_adapters```
 - Aggregate adapters: ```sh run_aggregate_adapters.sh -i /data/public/jkoubele/senescent_cells/detected_adapters -o /data/public/jkoubele/senescent_cells/aggregated_adapters```
   After that, check if the adapters were consistent across samples (the output will contain flag whether the adapters were detected successfully.)
 - Trimming reads: ```sh batch_trimming.sh -i /data/public/jkoubele/senescent_cells/FASTQ -o /data/public/jkoubele/senescent_cells/FASTQ_trimmed -a /data/public/jkoubele/cell_cultures_mtor/aggregated_adapters```
 - Alignment: ```sh batch_align.sh -i /data/public/jkoubele/senescent_cells/FASTQ_trimmed/ -o /data/public/jkoubele/senescent_cells/BAM -g /data/public/jkoubele/reference_genomes/GRCh38.p14```
 - Infer strandedness: ```sh run_infer_strandedness.sh -i /data/public/jkoubele/senescent_cells/BAM -o /data/public/jkoubele/senescent_cells/strandedness_info```
 - Feature counts (whole genes): ```sh batch_feature_counts.sh -i /data/public/jkoubele/senescent_cells/BAM -o /data/public/jkoubele/senescent_cells/feature_counts_genes -g /data/public/jkoubele/reference_genomes/GRCh38.p14 -a Homo_sapiens.GRCh38.112.gtf -s 2 -f gene```
 - Feature counts (exons only genes): ```sh batch_feature_counts.sh -i /data/public/jkoubele/senescent_cells/BAM -o /data/public/jkoubele/senescent_cells/feature_counts_exons -g /data/public/jkoubele/reference_genomes/GRCh38.p14 -a Homo_sapiens.GRCh38.112.gtf -s 2 -f exon```
 