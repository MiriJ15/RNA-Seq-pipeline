#!/bin/bash

# Path to the directory containing BAM files
BAM_DIR="alignments"

# Path to the GTF annotation file
ANNOTATION_FILE="GSE211993_gencode_fantomcat.v1.02.genes_transcripts_exons.gtf"

# Output directory for featureCounts results
OUTPUT_DIR="featureCounts_results"
mkdir -p "$OUTPUT_DIR"

# Output file name
OUTPUT_FILE="${OUTPUT_DIR}/featureCounts_results.txt"

# Run featureCounts
featureCounts -a "$ANNOTATION_FILE" \
              -o "$OUTPUT_FILE" \
              -F GTF \
              -t exon \
              -g gene_id \
              -p \
              --countReadPairs \
              --minOverlap 90 \
              --fracOverlap 0.9 \
              -T 32 \
              "$BAM_DIR"/*.bam
              
echo "featureCounts has completed."
