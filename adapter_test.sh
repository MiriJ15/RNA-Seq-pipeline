#!/bin/bash

# Adapter sequences
FORWARD_ADAPTER="AAGTCGGAGGCCAAGCGGTCTTAGGAAGACAA"
REVERSE_ADAPTER="AAGTCGGATCGTAGCCATGTCGTTCTGTGAGCCAAGGAGTTG"
FASTQ_DIR="path/to/fastq_files"

# Loop through all FASTQ files in the specified directory
for FASTQ_FILE in $FASTQ_DIR/*.fastq; do
    # Determine if the file is forward or reverse based on its name
    if [[ $FASTQ_FILE =~ _1.fastq$ ]]; then
        # Search for the forward adapter
        ADAPTER_SEQ=$FORWARD_ADAPTER
        COUNT=$(grep -o "$ADAPTER_SEQ" "$FASTQ_FILE" | wc -l)
        echo "Found $COUNT instances of the forward adapter sequence in $FASTQ_FILE"
    elif [[ $FASTQ_FILE =~ _2.fastq$ ]]; then
        # Search for the reverse adapter
        ADAPTER_SEQ=$REVERSE_ADAPTER
        COUNT=$(grep -o "$ADAPTER_SEQ" "$FASTQ_FILE" | wc -l)
        echo "Found $COUNT instances of the reverse adapter sequence in $FASTQ_FILE"
    else
        echo "Skipping $FASTQ_FILE: not recognized as forward or reverse."
    fi
done
