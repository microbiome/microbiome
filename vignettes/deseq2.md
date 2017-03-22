---
title: "DESeq2"
bibliography: 
- bibliography.bib
- references.bib
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
pseq <- dietswap
# Set baseline to 0 (in this data set it appears to be 1)
pseq <- microbiome::transform(pseq, "shift", shift = -1)
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
print(head(kable((df))))
```

```
## [1] "|     baseMean| log2FoldChange|     lfcSE|        stat|    pvalue|      padj|taxon                                 |"
## [2] "|------------:|--------------:|---------:|-----------:|---------:|---------:|:-------------------------------------|"
## [3] "| 1.179042e+02|     -5.6054637| 0.7137407|  -7.8536422| 0.0000000| 0.0000000|Uncultured Selenomonadaceae           |"
## [4] "| 6.462246e+03|     -3.7306691| 0.2855295| -13.0657932| 0.0000000| 0.0000000|Dialister                             |"
## [5] "| 1.006694e+03|     -3.4929604| 0.1878069| -18.5986829| 0.0000000| 0.0000000|Bacteroides intestinalis et rel.      |"
## [6] "| 8.769326e+04|     -2.9250398| 0.1987956| -14.7138053| 0.0000000| 0.0000000|Bacteroides vulgatus et rel.          |"
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
