#!/bin/bash

GENOME_FASTA="Homo_sapiens.GRCh38.dna.primary_assembl$
ANNOTATION_GTF="GSE211993_gencode_fantomcat.v1.02.gen$
STAR_INDEX_DIR="star_index"

mkdir -p "$STAR_INDEX_DIR"

# Run STAR in genomeGenerate mode
STAR --runThreadN 32 \
    --runMode genomeGenerate \
    --genomeDir "$STAR_INDEX_DIR" \
    --genomeFastaFiles "$GENOME_FASTA" \
    --sjdbGTFfile "$ANNOTATION_GTF" \  
    --sjdbOverhang 99 
    
echo "Indexing complete."
