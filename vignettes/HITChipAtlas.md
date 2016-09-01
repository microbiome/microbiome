---
title: "HITChip Atlas Overview"
author: "Leo Lahti, Jarkko Salojarvi, Anne Salonen, Willem M de Vos"
date: "2016-08-31"
bibliography: 
- bibliography.bib
- references.bib
output: 
  md_document:
    variant: markdown_github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial}
  %\usepackage[utf8]{inputenc}
-->




HITChip Atlas Overview
===========


The HITChip Atlas data collection covers microbial abundance profiles
for 130 genus-like taxa across over 10,000 samples, measured with 16S
rRNA phylogenetic microarrays. This page contains a reproducible
summary and overview of this data collection.

The analyses are based on the microbiome R package, which provides
convenient wrappers for many standard tasks encountered in microbiome
studies. The microbiome R package extends the
[phyloseq](http://joey711.github.io/phyloseq/import-data) class, a
standard representation format in R for taxonomic profiling.



### Data overview

Downloading the data in R:


```r
# Download the required R packages and then the HITChip Atlas data set
library("microbiome")

# Load the data (not yet public)
atlas <- readRDS("/home/lei/proj/hitchip-atlas-annotation-scripts/Atlas.RData")
```









