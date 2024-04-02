# ANGSD - RNAseq Project

## Overview
This project investigates the transcriptomic changes in neurons carrying PSEN1 mutations, focusing on the implications of these mutations in Alzheimer's disease pathology. The project leverages RNA sequencing data to explore changes in extracellular matrix components and non-coding gene expression associated with PSEN1 mutations.

## Data Sources
- **GEO Accession**: Data is available under GEO accession [GSE211993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211993).
- **SRA Study ID**: Raw sequence data can be found using ID [PRJNA873113](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA873113&o=acc_s%3Aa).
- **Publication**: The findings are detailed in the study titled "The transcriptomic landscape of neurons carrying PSEN1 mutations reveals changes in extracellular matrix components and non-coding gene expression," available [here](https://www.sciencedirect.com/science/article/pii/S0969996122003722?via%3Dihub).

## Methodology

### Prerequisites
- Installation of Conda and creation of a specific environment for the RNA-seq analysis is required. See `rna_seq_analysis.yaml` for environment setup.

### Data Preparation and Quality Control
- RNA extraction, library construction, and sequencing details are provided. Key tools for quality control include `fastqc`, `multiqc`, and `cutadapt`.

### Alignment and Feature Counting
- The project uses STAR for read alignment against the hg38 reference genome and featureCounts for gene expression quantification.

### Differential Expression Analysis
- Differential expression analysis focuses on comparing control samples against PSEN1 mutant samples, specifically analyzing the L150P mutation.

## Scripts and Analysis Pipeline
- Detailed scripts (found in `scripts` directory) for each analysis step are included, covering data downloading (`download_data.sh`), quality control (`fastqc.sh`, `multiqc.sh`), read alignment (`star.sh`), and differential expression analysis.
- A master script (`master.sh`) is provided to run all tasks sequentially for improved reproducibility.

## Results
- Results include differential expression analysis, PCA plots, volcano plots, and GO enrichment analysis, offering insights into the biological impact of PSEN1 mutations.
- Interactive plots and tables highlight significant genes and pathways affected by the mutations.

## Reproducibility
- The project ensures reproducibility through detailed documentation of the data sources, analysis pipeline, and the provision of scripts.

## Conclusion
TBC.
