- First QC: ```sh batch_qc.sh -i /data/public/jkoubele/celegans_mutants/FASTQ -o /data/public/jkoubele/celegans_mutants/QC_before_trimming```
- Detecting adapters: ```sh batch_detect_adapers.sh -i /data/public/jkoubele/celegans_mutants/FASTQ -o /data/public/jkoubele/celegans_mutants/detected_adapters```  
- No adapters were detected, also fastqc reported no adapter content. 
The reads in FASTQ differ in length and A content drops near read end -> the data seems to be already trimmed for adapters and poly A.  
- Alignment:  ```sh batch_align.sh -i /data/public/jkoubele/celegans_mutants/FASTQ -o /data/public/jkoubele/celegans_mutants/BAM -g /data/public/jkoubele/reference_genomes/WBcel235```
- Infer strandedness: ```sh run_infer_strandedness.sh -i /data/public/jkoubele/celegans_mutants/BAM -o /data/public/jkoubele/celegans_mutants/strandedness_info```
- Compute coverage: ```sh batch_compute_coverage.sh -i /data/public/jkoubele/celegans_mutants/BAM -o /data/public/jkoubele/celegans_mutants/coverage -g /data/public/jkoubele/reference_genomes/WBcel235 -f Caenorhabditis_elegans.WBcel235.dna.toplevel.fa.fai -s 2```
- Slopes estimation: ```sh batch_slope_estimation.sh -i /data/public/jkoubele/celegans_mutants/coverage -o /data/public/jkoubele/celegans_mutants/intron_slopes -g /data/public/jkoubele/reference_genomes/WBcel235```
- Intron selection: ```sh batch_select_introns.sh -i /data/public/jkoubele/celegans_mutants/intron_slopes -s /data/public/jkoubele/celegans_mutants/BAM -o /data/public/jkoubele/celegans_mutants/selected_intron_slopes```