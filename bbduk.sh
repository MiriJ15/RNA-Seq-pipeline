#!/bin/bash

# Directory paths
TRIMMED_FASTQ_DIR="trimmed_fastq"
RRNA_REMOVED_DIR="rRNA_removed_fastq"
TEMP_DIR="temp_rRNA_removed"

# Reference databases (ensure these are in FASTA format)
LSU_DB="SILVA_119_LSURef_tax_silva.fasta"
SSU_DB="SILVA_119.1_SSURef_Nr99_tax_silva.fasta"

# Number of threads to use
THREADS=32

mkdir -p "$RRNA_REMOVED_DIR" "$TEMP_DIR"

# Function to run BBDuk for rRNA removal
run_bbduk() {
    local in1=$1
    local in2=$2
    local out1=$3
    local out2=$4
    local ref=$5

    bbduk.sh \
        in="$in1" \
        in2="$in2" \
        out="$out1" \
        out2="$out2" \
        ref="$ref" \
        k=31 \
        ktrim=n \
        mcf=0.5 \
        tbo \
        tpe \
        threads="$THREADS"
}

# Loop through all the trimmed FASTQ files in pairs
for file1 in "$TRIMMED_FASTQ_DIR"/*_1_trimmed.fastq; do
    file2="${file1%_1_trimmed.fastq}_2_trimmed.fastq"
    
    # Intermediate and final output filenames
    temp1="$TEMP_DIR/$(basename "$file1")"
    temp2="$TEMP_DIR/$(basename "$file2")"
    final1="$RRNA_REMOVED_DIR/$(basename "$file1" .fastq)_rRNA_removed.fastq"
    final2="$RRNA_REMOVED_DIR/$(basename "$file2" .fastq)_rRNA_removed.fastq"
    
    # Step 1: Remove LSU rRNA
    echo "Removing LSU rRNA for: $(basename "$file1") and $(basename "$file2")"
    run_bbduk "$file1" "$file2" "$temp1" "$temp2" "$LSU_DB"
    
    # Step 2: Remove SSU rRNA using the output from step 1 as input
    echo "Removing SSU rRNA for: $(basename "$file1") and $(basename "$file2")"
    run_bbduk "$temp1" "$temp2" "$final1" "$final2" "$SSU_DB"
    
    # Optional: Remove intermediate files
    rm "$temp1" "$temp2"
    
    echo "SSU and LSU rRNA removal done for: $(basename "$file1") and $(basename "$file2")"
done

echo "rRNA removal process completed."
