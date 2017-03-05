---
title: "HITChip Atlas"
author: "Leo Lahti, Willem M de Vos et al."
date: "2017-03-05"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - atlas}
  %\usepackage[utf8]{inputenc}
-->


## Intestinal microbiota diversity in 1006 western adults

The data set from [Lahti et al. Nat. Comm. 5:4344,
2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) has microbiota profiling of 130 genus-like taxa across 1006 normal western adults from [Data Dryad](http://doi.org/10.5061/dryad.pk75d). Load the data in R:


```r
# Download the required R packages and then the HITChip Atlas data set
library(microbiome)
data(atlas1006)

# Let us give a generic name for this phyloseq data object,
# to be used further analysis
pseq <- atlas1006
```



### Estimating microbial diversity 


```r
library(phyloseq)
div <- microbiome::diversity(pseq, measures = c("Observed", "Shannon", "Simpson"))

library(knitr)
kable(head(div))
```



|         | Observed|  Shannon|   Simpson|
|:--------|--------:|--------:|---------:|
|Sample.1 |      130| 3.189726| 0.9230387|
|Sample.2 |      130| 3.396135| 0.9397719|
|Sample.3 |      130| 2.866104| 0.8850959|
|Sample.4 |      130| 3.058653| 0.9066459|
|Sample.5 |      130| 3.076850| 0.9184565|
|Sample.6 |      130| 2.945709| 0.8966565|


### Diversity vs. obesity (discrete factor)


```r
p <- plot_diversity(pseq, variable = "bmi_group", measures = c("Observed", "Shannon", "Simpson"), detection = 250, na.rm = TRUE)
print(p)
```

![plot of chunk div-example2](figure/div-example2-1.png)


