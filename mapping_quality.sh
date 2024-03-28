#!/bin/bash

# Define the directory containing STAR output files
STAR_OUTPUT_DIR="alignments"

# Define the directory where you want to save the MultiQC report
REPORT_DIR="$(pwd)/star_multiqc_report"

# Run MultiQC on the STAR outputs
echo "Running MultiQC on the STAR outputs..."
multiqc $STAR_OUTPUT_DIR -o $REPORT_DIR
