---
title: "HITChip Atlas Overview"
author: "Leo Lahti, Jarkko Salojarvi, Anne Salonen, Willem M de Vos"
date: "2016-06-15"
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
```

```
## Loading required package: phyloseq
```

```
## 
```

```
## 
## microbiome R package (microbiome.github.com)
##           
## 
## 
##  Copyright (C) 2011-2016
##           Leo Lahti and Jarkko Salojarvi 
## 
##         
##           <microbiome-admin@googlegroups.com>
```

```r
#data("atlas1006")
#pseq <- atlas1006
# Load the data (not yet public)
atlas = readRDS("/home/lei/proj/hitchip-atlas-annotation-scripts/Atlas.RData") # atlas
```

```
## Warning in gzfile(file, "rb"): cannot open compressed file '/home/lei/proj/
## hitchip-atlas-annotation-scripts/Atlas.RData', probable reason 'No such
## file or directory'
```

```
## Error in gzfile(file, "rb"): cannot open the connection
```









