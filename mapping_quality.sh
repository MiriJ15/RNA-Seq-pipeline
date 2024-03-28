#!/bin/bash

# Define directories
BAM_DIR="alignments"  # Directory containing BAM files
TEMP_DIR="flagstat_reports"  # Temporary directory for individual flagstat reports

# Create the temporary directory if it doesn't exist
mkdir -p "$TEMP_DIR"

# Initialize counters
total_reads=0
mapped_reads=0
properly_paired=0

# Loop through each BAM file in the directory
for BAM_FILE in "$BAM_DIR"/*.bam; do
    # Generate a temporary flagstat report
    TEMP_REPORT="${TEMP_DIR}/$(basename "$BAM_FILE").flagstat"
    samtools flagstat "$BAM_FILE" > "$TEMP_REPORT"
    
    # Extract and accumulate statistics
    while IFS= read -r line; do
        case "$line" in
            *in total*)
                total_reads=$(($total_reads + $(echo "$line" | cut -d ' ' -f 1)))
                ;;
            *mapped*)
                mapped_reads=$(($mapped_reads + $(echo "$line" | cut -d ' ' -f 1)))
                ;;
            *properly paired*)
                properly_paired=$(($properly_paired + $(echo "$line" | cut -d ' ' -f 1)))
                ;;
        esac
    done < "$TEMP_REPORT"
done

# Calculate percentages
mapped_percent=$(echo "scale=2; $mapped_reads / $total_reads * 100" | bc)
properly_paired_percent=$(echo "scale=2; $properly_paired / $total_reads * 100" | bc)

# Output combined statistics
echo "Combined statistics for all BAM files:"
echo "Total reads: $total_reads"
echo "Mapped reads: $mapped_reads ($mapped_percent%)"
echo "Properly paired reads: $properly_paired ($properly_paired_percent%)"

# Cleanup: Remove temporary directory and its contents
rm -r "$TEMP_DIR"
