---
title: "Example"
author: "Leo Lahti"
date: "2016-12-15"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - example}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


# Example document

This Rmarkdown document shows how to carry out some basic HITChip analysis with RStudio. [Download this file](https://raw.githubusercontent.com/microbiome/microbiome/master/vignettes/Example.Rmd), open it in RStudio, and then press the 'Knit HTML button'. This will generate and open a HTML file with some analysis with an example data set. You can then modify this file to use your own data (see below), or add new analyses.


## Update the microbiome package

Run this to make sure you have the latest version of the microbiome package:


```r
# Updating microbiome package
library(devtools)
install_github("microbiome/microbiome")
```

Load the tools and example data


```r
library(microbiome)
data("dietswap")
print(dietswap)
```

```
## phyloseq-class experiment-level object
## otu_table()   OTU Table:         [ 130 taxa and 222 samples ]
## sample_data() Sample Data:       [ 222 samples by 8 sample variables ]
## tax_table()   Taxonomy Table:    [ 130 taxa by 2 taxonomic ranks ]
```

