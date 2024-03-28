#!/bin/bash

# Define the number of threads to use
THREADS=32

# Path to the STAR genome index
GENOME_DIR="star_index"

# Directory containing rRNA removed FASTQ files
FASTQ_DIR="rRNA_removed_fastq"

# Output directory for the alignments
ALIGNMENT_DIR="alignments"
mkdir -p "$ALIGNMENT_DIR"

# Loop through the FASTQ files in pairs
for R1 in "$FASTQ_DIR"/*_1_trimmed_rRNA_removed.fastq; do
    # Infer the name of the R2 file from the name of the R1 file
    R2="${R1%_1_trimmed_rRNA_removed.fastq}_2_trimmed_rRNA_removed.fastq"

    # Extract the base name for output naming
    SAMPLE_NAME=$(basename "$R1" _1_trimmed_rRNA_removed.fastq)

    echo "Processing $SAMPLE_NAME"

    # Run STAR for alignment
    STAR --runMode alignReads \
         --runThreadN $THREADS \
         --genomeDir "$GENOME_DIR" \
         --readFilesIn "$R1" "$R2" \
         --outFileNamePrefix "${ALIGNMENT_DIR}/${SAMPLE_NAME}." \
         --outSAMtype BAM SortedByCoordinate

    echo "Alignment completed for $SAMPLE_NAME"
done

echo "All sample alignments have been completed."
