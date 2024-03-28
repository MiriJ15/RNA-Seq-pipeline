#!/bin/bash

# Adapter sequences
FORWARD_ADAPTER="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
REVERSE_ADAPTER="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"

# Create a directory for the trimmed output if it doesn't exist
mkdir -p trimmed_fastq

# Loop through all the _1.fastq files in the fastq directory
for file in fastq/*_1.fastq; do
    # Identify the pair
    file1=$file
    file2=${file/_1.fastq/_2.fastq}

    # Generate output filenames
    out1=trimmed_fastq/$(basename "$file1" .fastq)_trimmed.fastq
    out2=trimmed_fastq/$(basename "$file2" .fastq)_trimmed.fastq

    # Run Cutadapt
    cutadapt \
        -a "$FORWARD_ADAPTER" -A "$REVERSE_ADAPTER" \
        -q 30 --pair-filter=any --minimum-length=30 \
        -o "$out1" -p "$out2" \
        "$file1" "$file2"
done
