---
title: "Comparisons"
author: "Leo Lahti"
date: "2017-03-05"
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


## Group-wise comparisons

A number of methods for microbiota community comparisons have been proposed. For a recent benchmarking study, see [Weiss et al. (2017)](http://doi.org/10.1186/s40168-017-0237-y).

### Univariate comparisons

For individual taxa, diversity indicators etc.

 * [Linear mixed effect models](Mixedmodels.Rmd) 
 * [Negative binomial test](Negativebinomial.md) 

Other methods, not implemented here (see [Weiss et al. (2017)](http://microbiomejournal.biomedcentral.com/articles/10.1186/s40168-017-0237-y) for a recent survey):

 * [ANCOM](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4450248/); claimed to outperform ZIGs
 * [Zero-inflated Gaussians (ZIGs)](https://www.ncbi.nlm.nih.gov/pubmed/24076764/) (see [metagenomeSeq](https://bioconductor.org/packages/release/bioc/html/metagenomeSeq.html) Bioconductor package)
 * [DESeq2]() and other advanced methods based on negative binomial


### Multivariate comparisons

For community-level multivariate comparisons

 * [Multivariate linear models (limma)](limma.md)
 * [PERMANOVA](PERMANOVA.md)


