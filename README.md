# Karttunen & Patel et al. scripts

This is a repository of the source R code for the paper:

Karttunen, K., Patel, D., Xia, J. et al. Transposable elements as tissue-specific enhancers in cancers of endodermal lineage. Nature Communications 14, 5313 (2023). https://doi.org/10.1038/s41467-023-41081-4

#### This is the R code in Rmarkdown reports and associated utility scripts used to create the figures in the manuscript.
The scripts require the annotation files and preprocessed data as described in the methods.
All scripts have been tested on R 4.1.2 on Ubuntu 18.04.6 LTS.

#### All scripts require the following R packages:

GenomicRanges_1.44.0

ggrepel_0.9.2

ggsci 2.9

RColorBrewer 1.1-3

soGGi 1.24.1

tidyverse 1.3.2

wesanderson 0.3.6


## Figure 1.
Fig_1.rmd contains the R code to reproduce fig. 1 b, c, d and Supplementary figs. 1a, b, and c.
The script calculates the STARR-seq peak summit enrichment at TE classes lineages and creates the upset plot of the STARR peak and TE overlaps.

#### Required additional packages for Fig_1.rmd:

bedtoolsr 2.29.0-3

rstatix 0.7.1

UpSetR 1.4.0

## Figure 2.
Fig_2.rmd contains the R code to reproduce fig. 2 a, b c, d, Supplementary figs. 2a, b, c, 3c, and 4 a, b, c.
The script creates the enriched heatmap of the GP5d clusters, the motif enrichment heatmap, the TE enrichment per cluster and the GP5d WT vs GP5d p53-null peak enrichment.
The kmeans-clustering file required for the script is created with utils/kmeans-clustering.R.

#### Required additional packages for Fig_2.rmd:

bedtoolsr 2.29.0-3

BiocParallel 1.26.2

ChIPseeker 1.28.3

circlize 0.4.15

data.table 1.14.4

EnrichedHeatmap 1.22.0

ggpubr 0.4.0

matrixStats 0.62.0

rstatix 0.7.1

rtracklayer 1.52.1

TxDb.Hsapiens.UCSC.hg38.knownGene 3.13.0

Rsubread 2.6.4

## Figure 3.
Fig_3.rmd contains the R code to reproduce fig. 3 a, b, Supplementary figs. 3a, b, c, and 5a, b, c, d.
The script calculates the TE subfamily enrichment in GP5d and HepG2 and creates the subfamily-level motif heatmap.

#### Required additional packages for Fig_3.rmd:

bedtoolsr 2.29.0-3

circlize 0.4.15

ComplexHeatmap 2.8.0

ggpubr 0.4.0

rstatix 0.7.1

Rsubread 2.6.4

## Figure 5.
Fig_5.rmd contains the R code to reproduce figs. 5a, b, e, and f.
The script calculates the HepG2 non-methylated vs. methylated STARR-seq enrichment and creates the read count plots, the metaplot of signal enrichment and calculates the methylation differences between GP5d and HepG2.

#### Required additional packages for Fig_5.rmd:

bsseq 1.28.0

ggpubr 0.4.0

Rsubread 2.6.4

## Extended Data Fig. 5

This is a separate script for creating the figures in Extended Data Fig. 5. The cancer-specific TCGA ATAC peak data for the 23 cancer types were downloaded from https://gdc.cancer.gov/about-data/publications/ATACseq-AWG.

#### Required additional packages for Extended_data_fig_5.rmd:

circlize 0.4.15

ComplexHeatmap 2.8.0

## Figure 6.
Fig_6.rmd contains the R code to reproduce fig. 6a, b, and Supplementary figs. 9a, b.
The script processes the ABC and DEseq2 output data and creates the boxplots of TE contact gene expression.

RNA_seq_processing.rmd contains the scripts to process RNA-seq count data with DEseq2 to expression tables used in fig_6.Rmd.
Requires DESeq2 1.32.0 and tidyverse 1.3.2.
