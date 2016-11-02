---
title: "Density"
author: "Leo Lahti"
date: "2016-11-02"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - density}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->

### Abundance histograms

Population densities for Dialister:


```r
library(microbiome)
library(phyloseq)

# Example data from
# http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html
data(atlas1006)

# Pick the subset of RBB-preprocessed samples from time point 0
x <- subset_samples(atlas1006, time == 0 & DNA_extraction_method == "r")

# Visualize population densities for specific taxa
plot_density(x, "Dialister") + ggtitle("Absolute abundance")

# Same with log10 scaled X axis
plot_density(x, "Dialister", log10 = TRUE) + ggtitle("Log10")

# Same with log10 relative abundances
x <- transform_phyloseq(x, "relative.abundance")
tax <- "Dialister"
plot_density(x, tax, log10 = TRUE) + ggtitle("Relative abundance")
```

<img src="figure/hist-1.png" title="plot of chunk hist" alt="plot of chunk hist" width="280px" /><img src="figure/hist-2.png" title="plot of chunk hist" alt="plot of chunk hist" width="280px" /><img src="figure/hist-3.png" title="plot of chunk hist" alt="plot of chunk hist" width="280px" />

