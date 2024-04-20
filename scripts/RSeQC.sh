#!/bin/bash

# Define the directory containing BAM files and where to output results
BAM_DIR="alignments"
OUTPUT_DIR="rRNA_overlap_results"
RIBO_BED="rRNA.bed"

# Create output directory if it does not exist
mkdir -p "$OUTPUT_DIR"

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Extract the base name of the file for output naming
    BASE_NAME=$(basename "$BAM_FILE" .bam)
    
    # Perform bedtools intersect
    bedtools intersect -abam "$BAM_FILE" -b "$RIBO_BED" -bed -wa | wc -l > "$OUTPUT_DIR/${BASE_NAME}_rRNA_count.txt"
    
    echo "Processed $BAM_FILE"
done

echo "All files processed. rRNA read counts are in $OUTPUT_DIR."