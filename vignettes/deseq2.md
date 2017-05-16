---
title: "DESeq2"
bibliography: 
- bibliography.bib
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Normalization and group-wise comparisons with DESeq2

Examples adapted from [Callahan et al. F1000 (2017)](https://f1000research.com/articles/5-1492/v2).

Load example data:


```r
# Load libraries
library(microbiome)
library(ggplot2)

# Probiotics intervention example data 
data(dietswap) 

# Set baseline to 0 (in this data set it appears to be 1)
pseq <- microbiome::transform(dietswap, "shift", shift = -1)

# Only check the core taxa to speed up examples
pseq <- core(pseq, detection = 10^3, prevalence = 95/100)
```


Toy example, to be polished:


```r
library(phyloseq)
library(structSSI)
library(plyr)
library(dplyr)
library(reshape2)
library(DESeq2)

# Running the DESeq2 analysis
ds2 <- phyloseq_to_deseq2(pseq, ~ nationality)
dds <- DESeq(ds2)
res <- results(dds)
df <- as.data.frame(res)
df$taxon <- rownames(df)
df <- df %>% arrange(log2FoldChange, padj)

library(knitr)
print(head(kable((df))))
```

```
## [1] "|   baseMean| log2FoldChange|     lfcSE|        stat|    pvalue|      padj|taxon                                |"
## [2] "|----------:|--------------:|---------:|-----------:|---------:|---------:|:------------------------------------|"
## [3] "|  88896.518|     -3.0629202| 0.1999330| -15.3197368| 0.0000000| 0.0000000|Bacteroides vulgatus et rel.         |"
## [4] "|  16380.395|     -2.6968824| 0.1546172| -17.4423188| 0.0000000| 0.0000000|Allistipes et rel.                   |"
## [5] "|   5466.663|     -1.4142640| 0.1443252|  -9.7991478| 0.0000000| 0.0000000|Bryantella formatexigens et rel.     |"
## [6] "|  13000.276|     -0.5528212| 0.1291264|  -4.2812424| 0.0000186| 0.0000351|Subdoligranulum variable at rel.     |"
```


Validating DESeq2 results


```r
# Identify top taxa based on standard ANOVA
source(system.file("extdata/check_anova.R", package = "microbiome"))
ano <- check_anova(pseq, "nationality");
ano$log2FC <- log2(ano$ave.AFR) - log2(ano$ave.AAM)
taxa.anova <- as.character(subset(ano, padj < 0.01 & abs(log2FC) > log2(2))$taxa)

# Pick the top taxa based on DESEq2
taxa.deseq <- subset(res.deseq, padj < 0.01 & abs(log2FoldChange) > log2(2))$taxon

# Check overlap
# Most DESEq2 taxa are confirmed with ANOVA
library(gplots)
venn( list(ANOVA = taxa.anova,DESeq2 = taxa.deseq) )

# Also the lowest p-values are well correlated (higher not so)
plot(log10(res.deseq$padj), log10(ano$padj), xlab = "DESeq2 adjusted p-value", ylab("ANOVA adjusted p-value"))
```
