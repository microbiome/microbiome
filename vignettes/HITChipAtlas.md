---
title: "HITChip Atlas Overview"
author: "Leo Lahti, Jarkko Salojarvi, Anne Salonen, Willem M de Vos"
date: "2016-06-08"
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

We use the [phyloseq](http://joey711.github.io/phyloseq/import-data)
class, a standard representation format in R for taxonomic
profiling. This package provides extensions and convenient wrappers
for many standard tasks encountered in microbiome studies. 


## Intestinal microbiota diversity 

This data set has microbiota profiling of 130 genus-like taxa across over 10,000 samples.

Download the data in R:


```r
# Download the required R packages and then the HITChip Atlas data set
library("microbiome")
#data("atlas1006")
#pseq <- atlas1006
# Load the data (not yet public)
atlas = readRDS("/home/lei/proj/hitchip-atlas-annotation-scripts/Atlas.RData") # atlas
```


### Diversity vs. age


```r
# Pick the subset of RBB-preprocessed samples from time point 0
pseq <- subset_samples(atlas, time == 0 & DNA_extraction_method == "rbb")

# Visualize
p <- plot_regression(diversity~age, sample_data(pseq))
print(p)
```

![plot of chunk atlas-example3](figure/atlas-example3-1.png)



### Licensing and Citations

This material can be freely used, modified and distributed under the
[Two-clause FreeBSD
license](http://en.wikipedia.org/wiki/BSD\_licenses). Kindly cite as
'Leo Lahti, Jarkko Salojarvi, Anne Salonen and Willem M de
Vos. HITChip Atlas. URL: http://microbiome.github.com'.


### References





