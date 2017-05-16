---
title: "Community comparisons"
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


## Group-wise comparisons

A number of methods for microbiota community comparisons have been proposed. For a recent benchmarking study, see [Weiss et al. (2017)](http://doi.org/10.1186/s40168-017-0237-y). For a comprehensive example workflow, see [Callahan et al. F1000 (2017)](https://f1000research.com/articles/5-1492/v2).

### Univariate comparisons

For individual taxa, diversity indicators etc.

 * [Linear mixed effect models](Mixedmodels.html) 
 * [Negative binomial test](Negativebinomial.html)
 * [ANCOM](ANCOM.html)  

Other methods, not implemented here (see [Weiss et al. (2017)](http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) for a recent survey):

 * [Zero-inflated Gaussians (ZIGs)](https://www.ncbi.nlm.nih.gov/pubmed/24076764/) (see [metagenomeSeq](https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html) Bioconductor package)
 * [DESeq2](deseq2.Rmd) and other advanced methods based on negative binomial


### Multivariate comparisons

For community-level multivariate comparisons

 * [Multivariate linear models (limma)](limma.html)
 * [PERMANOVA](PERMANOVA.html)

