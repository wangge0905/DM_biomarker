# DM_biomarker

R scripts for plasma EV-associated RNA analysis in anti-MDA5-positive and anti-ARS-positive IIM-ILD, including RNA profiling, subtype classification and pulmonary-impairment severity analysis.

## Contents

- `inputs/`: RNA count matrices, sample metadata, gene annotations and Hallmark gene sets.
- `folds/`: training-set cross-validation folds and severity-model inner and outer folds.
- `results/`: selected RNA panels, participant-level predictions and performance summaries.

## Scripts

- **RNA profiling:** `QC.R`, `QC_region.R`, `RNA_biotype.R`, `cf-biotype.R`, `cf-read-ratio.R`, `TE.R`, `TE_heatmap.R`.
- **Differential expression:** `ComBat_seq.R`, `differential_expression.R`, `cf-differential.R`, `cf-differential-plot.R`.
- **Subtype classification:** `biomarker.R`, `biomarker_panel.R`, `biomarker_permutation.R`, `selected_gene_expression.R`.
- **Severity analysis:** `ev-pathway-score.R`, `severity_RNA.R`, `severity_pathway.R`, `severity_ROC.R`, `severity_plot.R`, `severity_confusion.R`, `severity_clinical_plot.R`.
- **Sensitivity and permutation analyses:** `age_sex_sensitivity.R`, `clinical_sensitivity.R`, `clinical_sensitivity_plot.R`, `severity_rank.R`, `severity_permutation.R`, `severity_permutation_plot.R`.
- **miRNA tissue expression:** `miRNA_organ.R`.

## Software

R 4.3.2. Package versions used:

```text
data.table    1.15.2
dplyr         1.1.4
tidyr         1.3.1
stringr       1.5.1
reshape2      1.4.4
openxlsx      4.2.8
edgeR         4.0.16
sva           3.35.2
e1071         1.7-14
glmnet        4.1-9
ggplot2       3.5.0
pheatmap      1.0.12
ggpubr        0.6.0
ggrepel       0.9.5
ggsci         3.0.1
ggsignif      0.6.4
patchwork     1.2.0
cowplot       1.1.3
gtable        0.3.4
scales        1.3.0
RColorBrewer  1.1-3
ragg          1.2.7
```
