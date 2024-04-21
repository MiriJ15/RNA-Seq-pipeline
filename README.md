# ANGSD - RNAseq Project

## Overview
This project explores transcriptomic alterations in neurons with PSEN1 mutations, to understand their role in Alzheimer's disease pathology. Using RNA sequencing data, we analyze changes in extracellular matrix components and non-coding gene expression linked to these mutations.

## Data Sources
- **GEO Accession**: Data available under [GSE211993](https://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE211993).
- **SRA Study ID**: Access raw sequence data with ID [PRJNA873113](https://www.ncbi.nlm.nih.gov/Traces/study/?acc=PRJNA873113&o=acc_s%3Aa).
- **Publication**: Find detailed study insights [here](https://www.sciencedirect.com/science/article/pii/S0969996122003722?via%3Dihub).

## Methodology

### Prerequisites
- Requires Conda environment setup. Refer to `rna_seq_analysis.yaml` for environment details.

### Data Preparation and Quality Control
- Procedures include RNA extraction, library construction, and sequencing. Key QC tools: `fastqc`, `multiqc`, `cutadapt`.

### Alignment and Feature Counting
- Utilizes STAR for alignment against the hg38 genome and featureCounts for expression quantification.

### Differential Expression Analysis
- Focuses on comparing control and PSEN1 mutant samples, particularly analyzing the L150P mutation.

## Scripts and Analysis Pipeline
- Scripts for each analysis step are in the `scripts` directory, including data download (`download_data.sh`), QC (`fastqc.sh`, `multiqc.sh`), alignment (`star.sh`), and differential expression analysis (`project.Rmd`/`project.html`).
- `master.sh` runs all scripts sequentially for enhanced reproducibility.

## Results
- Includes differential expression analysis, PCA plots, volcano plots, and GO enrichment analysis, highlighting the biological impact of PSEN1 mutations.
- Interactive plots and tables emphasize significant genes and pathways.

## Reproducibility
- Ensures reproducibility through detailed documentation of data sources, analysis pipeline, and scripts.

## Conclusion
- The findings reinforce the crucial role of PSEN1 mutations in Alzheimer's disease, providing insights that could lead to novel therapeutic targets. Final conclusions will detail the implications of these transcriptomic changes in disease progression.

---

## Additional Project Details

### Introduction
This project addresses how mutations in the PSEN1 gene affect gene expression in neurons and contribute to Alzheimer's disease pathology. It focuses on the L150P mutation known for altering amyloid precursor protein processing, which is pivotal in the disease's early onset.

### Hypothesis
We hypothesize that mutations in the PSEN1 gene result in significant changes in the expression of genes associated with the extracellular matrix, influencing Alzheimer's disease pathogenesis. This analysis aims to identify these critical gene expression alterations using RNA-seq data.

### Key Insights from RNA-seq Analysis
Initial findings from RNA-seq have identified numerous genes significantly altered by the PSEN1 mutation, with a focus on those involved in extracellular matrix organization. This supports the hypothesis that changes in the extracellular matrix play a significant role in the disease's progression.

### Future Directions
Future analyses will expand on these findings by incorporating a broader array of samples and employing techniques like single-cell RNA-seq to delve deeper into the cellular impacts of PSEN1 mutations.

### Limitations
The current study's limitations include the homogeneous nature of the sample set and the potential batch effects due to processing. Future studies will aim to mitigate these limitations by diversifying the sample pool and refining experimental protocols.

## Citations
[1] Corsi, G. I., Gadekar, V. P., Haukedal, H., Doncheva, N. T., Anthon, C., Ambardar, S., Palakodeti, D., Hyttel, P., Freude, K., Seemann, S. E., & Gorodkin, J. (2023, March). The transcriptomic landscape of neurons carrying PSEN1 mutations reveals changes in extracellular matrix components and non-coding gene expression. Neurobiology of Disease, 178, 105980. https://doi.org/10.1016/j.nbd.2022.105980

[2] Bagaria, J., Bagyinszky, E., & An, S. S. A. (2022, September 19). Genetics, Functions, and Clinical Impact of Presenilin-1 (PSEN1) Gene. International Journal of Molecular Sciences, 23(18), 10970. https://doi.org/10.3390/ijms231810970

[3] Poon, A., Schmid, B., Pires, C., Nielsen, T. T., Hjermind, L. E., Nielsen, J. E., Holst, B., Hyttel, P., & Freude, K. K. (2016, November). Generation of a gene-corrected isogenic control hiPSC line derived from a familial Alzheimer’s disease patient carrying a L150P mutation in presenilin 1. Stem Cell Research, 17(3), 466–469. https://doi.org/10.1016/j.scr.2016.09.018

## Key Data Sets
| File                                     | Location        | Description                                                       |
|------------------------------------------|-----------------|-------------------------------------------------------------------|
| counts_for_each_replicate.png            | key_data_sets   | Image showing the count for each replicate, categorized by status.|
| biotype_counts.txt                       | key_data_sets   | Text file listing the number of each gene biotype in the significant genes dataset. |
| top_differentially_expressed_genes.csv   | key_data_sets   | CSV file detailing prevalent protein-coding gene IDs and their gene names. |
| significant_fold_change_genes.csv        | key_data_sets   | CSV file documenting all significant genes with substantial fold change. |
| enriched_BP_GO_terms.csv                 | key_data_sets   | CSV file containing Biological Process GO terms from the dataset. |
| enriched_CC_GO_terms.csv                 | key_data_sets   | CSV file listing Cellular Component GO terms identified. |
| enriched_MF_GO_terms.csv                 | key_data_sets   | CSV file detailing Molecular Function GO terms discovered. |
