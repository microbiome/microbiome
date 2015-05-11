---
title: "Atlas workflow"
author: "Leo Lahti"
date: "2015-05-11"
output:
  md_document:
    toc: false
    variant: markdown_github
  html_document:
    toc: false
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{Project Template}
  %\usepackage[utf8]{inputenc}
-->


## Intestinal microbiota diversity in 1006 western adults

We investigate here an example data set from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) that contains large-scale profiling of 130 genus-like taxa across 1006 normal western subjects. This data set is available for download from the open [Data Dryad](http://doi.org/10.5061/dryad.pk75d) repository.


### Download HITChip Atlas data

[Load the HITChip Atlas microbiome profiling data in R](Data.md)


```r
# Download the required R packages and then the HITChip Atlas data set
library("rdryad")
library("microbiome")
```

```
## Loading required package: ade4
## Loading required package: dplyr
## 
## Attaching package: 'dplyr'
## 
## The following object is masked from 'package:stats':
## 
##     filter
## 
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
## 
## Loading required package: e1071
## Loading required package: phyloseq
## Loading required package: reshape2
## Loading required package: vegan
## Loading required package: permute
## Loading required package: lattice
## This is vegan 2.2-1
## 
## Attaching package: 'vegan'
## 
## The following object is masked from 'package:ade4':
## 
##     cca
```

```
## Warning: replacing previous import by 'dplyr::select' when loading
## 'microbiome'
```

```
## 
## microbiome R package (microbiome.github.com)
##           
## 
## 
##  Copyright (C) 2011-2015
##           Leo Lahti and Jarkko Salojarvi 
## 
##         
##           <microbiome-admin@googlegroups.com>
## 
## 
## Attaching package: 'microbiome'
## 
## The following object is masked from 'package:lattice':
## 
##     densityplot
## 
## The following object is masked from 'package:e1071':
## 
##     impute
```

```r
d <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```


### Diversity 


```r
# Estimat diversity using the vegan package
# NOTE: data needs to be in absolute scale, not logarithmic!
# Also pay attention on having samples on the rows (not on columns)
di <- vegan::diversity(atlas$microbes, index = "shannon")
```

```
## Error in as.matrix(x): object 'atlas' not found
```

```r
hist(di, main = "Diversity")
```

```
## Error in hist(di, main = "Diversity"): object 'di' not found
```

### Compare with known background factors

Diversity vs. obesity


```r
par(mar = c(6, 4, 3, 1))
bmi <- atlas$meta$BMI_group
```

```
## Error in eval(expr, envir, enclos): object 'atlas' not found
```

```r
di <- vegan::diversity(atlas$microbes)
```

```
## Error in as.matrix(x): object 'atlas' not found
```

```r
boxplot(di ~ bmi, las = 2, main = "Microbiota diversity vs. obesity")
```

```
## Error in eval(expr, envir, enclos): object 'di' not found
```


Diversity vs. age with smoothed confidence intervals:


```r
library(microbiome)
library(sorvi)
library(dplyr)

# Pick the subset of RBB-preprocessed samples from time point 0
rbb.samples <- filter(atlas$meta, Time == 0 & DNA_extraction_method == "r")$SampleID
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'atlas' not found
```

```r
# Collect variables into a data frame
df <- data.frame(Age = atlas$meta[rbb.samples, "Age"], Diversity = di[rbb.samples])
```

```
## Error in data.frame(Age = atlas$meta[rbb.samples, "Age"], Diversity = di[rbb.samples]): object 'atlas' not found
```

```r
# Visualize
p <- sorvi::regression_plot(Diversity~Age, df, shade = TRUE, mweight = TRUE, verbose = FALSE)
```

```
## Error in data[[IV]]: object of type 'closure' is not subsettable
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```

