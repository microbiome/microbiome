---
title: "Linear models"
author: "Leo Lahti"
date: "2016-11-13"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->

## Linear model examples

### Continuous variables

Rapid quantification of continuous associations can be done with the
lm_phyloseq wrapper function. This uses limma model to generate a
table of P-values and effect sizes. No confounding variables taken
into account in this wrapper. See the [limma
homepage](http://bioinf.wehi.edu.au/limma/) for more info.


```r
library(limma)
library(microbiome)
data("atlas1006")
# Pick RBB extracted samples (r) and baseline time point
pseq <- subset_samples(atlas1006, DNA_extraction_method == "r" & time == 0)
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
tab <- lm_phyloseq(pseq, "age")
kable(head(tab), digits = 3)
```



|                                 |  logFC| AveExpr|      t| P.Value| adj.P.Val|      B|
|:--------------------------------|------:|-------:|------:|-------:|---------:|------:|
|Clostridium orbiscindens et rel. |  0.005|   4.122|  4.542|       0|     0.001|  0.696|
|Eggerthella lenta et rel.        |  0.003|   2.596|  4.451|       0|     0.001|  0.306|
|Butyrivibrio crossotus et rel.   |  0.004|   4.039|  4.232|       0|     0.001| -0.604|
|Bifidobacterium                  | -0.007|   4.034| -4.166|       0|     0.001| -0.872|
|Peptococcus niger et rel.        |  0.003|   2.294|  3.872|       0|     0.003| -2.010|
|Eubacterium hallii et rel.       |  0.004|   3.661|  3.688|       0|     0.005| -2.682|


