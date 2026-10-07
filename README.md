# Integrated expression, mutation, and survival analysis of 17 key genes in breast cancer using TCGA-BRCA data

This repository contains the R scripts used for the analyses and figures reported in our published research article:

**Integrated expression, mutation, and survival analysis of 17 key genes in breast cancer using TCGA-BRCA data**

Sophia Wang, Christopher Wang, and Haiyan Lei  
*Journal of Emerging Investigators*, 2026  
Published: August 6, 2026  
DOI: https://doi.org/10.59720/25-295

## Overview

This study used publicly available data from The Cancer Genome Atlas Breast Cancer cohort (TCGA-BRCA) to investigate gene expression, somatic mutations, and patient survival for 17 genes with established relevance to breast cancer.

The genes analyzed were:

- **Tumor suppressor genes:** BRCA1, BRCA2, TP53, PTEN, RB1, CHEK2, CDH1, MAP3K1
- **Oncogenes:** PIK3CA, AKT1, ERBB2, MYC, CCND1, FGFR1
- **Hormone receptor-related genes:** ESR1, PGR
- **Transcription factor:** GATA3

The analyses included:

- RNA-seq preprocessing and differential expression analysis
- Volcano plot visualization
- Gene-expression heatmap analysis
- Kaplan-Meier survival analysis
- Cox proportional hazards regression
- Somatic mutation analysis
- Comparison of gene expression between wild-type and mutant tumors

## Repository Contents

| File | Description |
| --- | --- |
| `01_differential_expression.R` | TCGA-BRCA RNA-seq preprocessing and differential expression analysis |
| `02_figure2_volcano_plot.R` | Generates Figure 2: differential-expression volcano plot |
| `03_figure3_heatmap.R` | Generates Figure 3: heatmap of expression patterns for the 17 selected genes |
| `04_figure4_survival_analysis.R` | Generates Figure 4: Kaplan-Meier survival curves and Cox proportional hazards analyses |
| `05_figure5_mutation_analysis.R` | Generates Figure 5: mutation landscape and wild-type vs. mutant expression analysis for TP53, BRCA1, and BRCA2 |

## Data

All data analyzed in this study are publicly available through The Cancer Genome Atlas (TCGA) via the Genomic Data Commons (GDC).

The TCGA-BRCA cohort was used for RNA-seq gene expression, somatic mutation, and clinical survival data.

The analyses reported in the paper included:

- **Expression analysis:** 991 tumor samples and 113 normal samples after preprocessing and quality control
- **Survival analysis:** 850 tumor samples with clinical follow-up data
- **Mutation analysis:** 878 tumor samples with available somatic mutation data

Because TCGA data availability differs among data types, not every patient/sample was included in every analysis.

## Analysis

### Differential Expression

Raw RNA-seq count data were retrieved using `TCGAbiolinks`.

Gene-level counts were normalized and analyzed using `DESeq2`. Low-count genes were filtered before differential expression analysis.

Differential expression between tumor and normal samples was evaluated using the DESeq2 Wald test.

Genes with:

- absolute log2 fold change > 1
- adjusted p-value < 0.05

were considered statistically significant.

### Heatmap

Expression patterns for the 17 selected genes were visualized using row-wise Z-score normalization.

For visualization, Z-scores were clipped to the range -3 to +3.

### Survival Analysis

Overall survival was calculated using TCGA clinical data.

Patients were divided into high- and low-expression groups according to median gene expression.

Survival analyses included:

- Kaplan-Meier survival curves
- log-rank tests
- Cox proportional hazards regression
- hazard ratios with 95% confidence intervals

Kaplan-Meier curves were generated for BRCA1, BRCA2, and TP53, while Cox regression was performed for all 17 genes.

### Mutation Analysis

Somatic mutation data in Mutation Annotation Format (MAF) were analyzed using `maftools`.

Mutation analyses focused on:

- TP53
- BRCA1
- BRCA2

RNA-seq-derived transcript abundance was compared between wild-type and mutant tumors using the Wilcoxon rank-sum test.

## Software and R Packages

The analyses reported in the publication were performed using **R version 4.5.1**.

Major packages included:

- `TCGAbiolinks`
- `DESeq2`
- `maftools`
- `survival`
- `survminer`
- `ggplot2`
- `EnhancedVolcano`
- `pheatmap`

Additional package dependencies may be loaded by the individual scripts.

## Key Findings

The study identified substantial transcriptional differences between breast tumor and normal tissues.

Among the findings:

- ERBB2, ESR1, and CCND1 showed increased expression in tumors.
- BRCA2 and CHEK2 showed increased transcript abundance in tumors.
- MYC and RB1 showed reduced mRNA expression.
- TP53 was the most frequently mutated of the three genes examined in the mutation analysis.
- Expression levels of BRCA1, BRCA2, and TP53 showed trends in survival analyses, but these associations did not reach statistical significance.

These findings emphasize that gene expression, mutation status, and clinical outcomes are not necessarily directly aligned in breast cancer.

## Reproducibility

The scripts in this repository correspond to the final analyses used in the published article.

The scripts are numbered according to the general analysis workflow and figure order:

1. Differential expression analysis
2. Figure 2 — volcano plot
3. Figure 3 — expression heatmap
4. Figure 4 — survival analysis
5. Figure 5 — mutation analysis

TCGA data are not redistributed in this repository. The scripts obtain and/or process publicly available TCGA-BRCA data as described in the article and individual scripts.

## Citation

If you use or reference this analysis, please cite:

> Wang S, Wang C, Lei H. Integrated expression, mutation, and survival analysis of 17 key genes in breast cancer using TCGA-BRCA data. Journal of Emerging Investigators. 2026. https://doi.org/10.59720/25-295

## Article

The published article is available at:

https://doi.org/10.59720/25-295

## Authors

- Sophia Wang
- Christopher Wang
- Haiyan Lei

## License

The published article is distributed under the license specified by the Journal of Emerging Investigators.

The source code in this repository is provided for research and educational reproducibility. No separate software license has currently been specified for the repository.
