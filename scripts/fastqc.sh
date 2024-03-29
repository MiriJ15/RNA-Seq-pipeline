#!/bin/bash

FASTQ_DIR="fastq"

# Output directory for FastQC results
OUTPUT_DIR="fastqc_results"

# Create output directory if it doesn't exist
mkdir -p "$OUTPUT_DIR" multiqc

# Run FastQC on all FASTQ files in the directory
fastqc -o "$OUTPUT_DIR" -t 40 "$FASTQ_DIR"/*.fastq

# Run MultiQC to combine the FastQC reports
multiqc "$OUTPUT_DIR" -o multiqc
