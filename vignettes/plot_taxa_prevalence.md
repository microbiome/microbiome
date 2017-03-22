---
title: "Plot Taxa Prevalence"
author: "Sudarshan A. Shetty"
date: "2017-03-05"
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - Plot Taxa Prevalence}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Load data

Load [example data](Data.md):  
For this example we will use data from [Halfvarson J., et al. Nature Microbiology, 2017](http://www.nature.com/articles/nmicrobiol20174). It was downloaded from [Qitta](https://qiita.ucsd.edu/study/description/1629).  



```r
library(microbiome)
data(DynamicsIBD)
p0 <- DynamicsIBD
```

## Plot Taxa Prevalence

Phylum level:


```r
p0.f <- format_phyloseq(p0)
plot <- plot_taxa_prevalence(p0.f, 'Phylum')
print(plot)
```

![plot of chunk unnamed-chunk-2](figure/unnamed-chunk-2-1.png)

You can also plot the prevalence at Order/Class/Family level.  

