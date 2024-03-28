#!/bin/bash

# Downloading the SILVA LSU and SSU reference FASTA files
wget https://ftp.arb-silva.de/release_119.1/Exports/SILVA_119_LSURef_tax_silva.fasta.gz
wget https://ftp.arb-silva.de/release_119.1/Exports/SILVA_119.1_SSURef_Nr99_tax_silva.fasta.gz

# Decompressing the downloaded files
gunzip SILVA_119_LSURef_tax_silva.fasta.gz
gunzip SILVA_119.1_SSURef_Nr99_tax_silva.fasta.gz
