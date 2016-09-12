---
title: "HITChip Atlas Overview"
author: "Leo Lahti, Jarkko Salojarvi, Anne Salonen, Willem M de Vos"
date: "2016-09-13"
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

The data contains:

 * ``10763`` samples
 * ``5030`` unique subjects
 * ``51.8``% female ratio. Altogether, sex information is available for ```round(100 * sum(table(subset(sample_data(atlas))$sex)[c("male", "female")])/nsamples(atlas), 1)```% of the samples.
 * ``3810`` samples (``35.4``%) from ``1798`` unique subjects with reported health problems ('compromised')
 * ``455`` samples from ``208`` unique subjects with reported antibiotics use.
 * ``112`` samples from ``77`` subjects with reported probiotics use.
 * ``1150`` samples from ``813`` subjects with reported medication.  
 

<img src="figure/hatlas-sampletype-1.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="280px" /><img src="figure/hatlas-sampletype-2.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="280px" /><img src="figure/hatlas-sampletype-3.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="280px" /><img src="figure/hatlas-sampletype-4.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="280px" /><img src="figure/hatlas-sampletype-5.png" title="plot of chunk hatlas-sampletype" alt="plot of chunk hatlas-sampletype" width="280px" />


### Diversity vs. host variables


```r
# Pick the subset of RBB-preprocessed samples from time point 0
pseq <- subset_samples(atlas, time == 0 & DNA_extraction_method == "rbb")

p <- plot_regression(diversity~age, sample_data(pseq))
print(p)
```

<img src="figure/hatlas-example3-1.png" title="plot of chunk hatlas-example3" alt="plot of chunk hatlas-example3" width="800px" />

```r
#p <- boxplot_abundance(pseq, x = "bmi_group", y = "Dialister")
#print(p)
```



### Licensing and Citations

This material can be freely used, modified and distributed under the
[Two-clause FreeBSD
license](http://en.wikipedia.org/wiki/BSD\_licenses). Kindly cite as
'Leo Lahti, Jarkko Salojarvi, Anne Salonen and Willem M de
Vos. HITChip Atlas. URL: http://microbiome.github.com'.


### References





