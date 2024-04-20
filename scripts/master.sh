#!/bin/bash

# Stop on error
set -e

# Activate environment
conda activate rna_seq_analysis

# Step 1: Download Data
bash download_data.sh

# Step 2: Quality Control
bash fastqc.sh

# Step 3: Adapter Trimming
bash cutadapt.sh

# Step 4: rRNA Removal
bash silva_database_download.sh
bash run_bbduk.sh

# Quality control after processing
bash processed_multiqc.sh

# Step 5: Alignment
bash hg38.sh
bash index.sh
bash star.sh

# Quality check for mapping
bash mapping_quality.sh

# Assess rRNA removal via RSeQC
bash RSeQC.sh

# Step 6: Feature Counts
bash featurecounts.sh

echo "Analysis Complete"
