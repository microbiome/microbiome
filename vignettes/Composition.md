---
title: "Microbiome composition"
bibliography: 
- bibliography.bib
- references.bib
output: 
  prettydoc::html_pretty:
    theme: cayman
    highlight: github
---
  <!--
  %\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{microbiome tutorial - composition}
%\usepackage[utf8]{inputenc}
%\VignetteEncoding{UTF-8}  
-->
  
  
## Microbiota composition
  

Also see [phyloseq barplot examples](http://joey711.github.io/phyloseq/plot_bar-examples.html).
  
Read example data from a [diet swap study](http://dx.doi.org/10.1038/ncomms7342):
  

```r
# Example data
library(microbiome)
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
##  Copyright (C) 2011-2017 Leo Lahti et al. <microbiome.github.io>
```

```
## 
## Attaching package: 'microbiome'
```

```
## The following object is masked from 'package:base':
## 
##     transform
```

```r
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
data(dietswap)

# Just use prevalent taxa to speed up examples
pseq <- core(dietswap, detection = 0.1/100, prevalence = 90)
```

```
## Error in validObject(.Object): invalid class "otu_table" object: 
##  OTU abundance data must have non-zero dimensions.
```

```r
# Pick sample subset
library(phyloseq)
pseq2 <- subset_samples(pseq, group == "DI" & nationality == "AFR")
```

```
## Error in sample_data(physeq): object 'pseq' not found
```

### Barplots for composition

Show OTU absolute abundance in each sample. Plot absolute taxon
abundances:
  

```r
library(ggplot2)
theme_set(theme_bw(22)) # Black/white color theme
p <- plot_composition(pseq2, taxonomic.level = "Phylum") +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 10, byrow = TRUE))
```

```
## Error in aggregate_taxa(x, taxonomic.level): object 'pseq2' not found
```

```r
print(p)       
```

```
## Error in print(p): object 'p' not found
```

Arrange by sample variable. Focus on the core taxa. Africans have more Prevotella as expected. Absolute counts:
  

```r
# Limit the analysis on core taxa and specific sample group
pseq2 <- pseq %>%
  core(detection = 10^4, prevalence = .5) %>%
  subset_samples(group == "DI" & timepoint.within.group == 1)
```

```
## Error in eval(expr, envir, enclos): object 'pseq' not found
```

```r
p <- plot_composition(pseq2,
                      sample.sort = "nationality", # Sort by nationality
                      x.label = "nationality") +   # Label by nationality
  guides(fill = guide_legend(nrow = 5, byrow = TRUE)) +
  theme(legend.position = "bottom")
```

```
## Error in plot_composition(pseq2, sample.sort = "nationality", x.label = "nationality"): object 'pseq2' not found
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```


Same with compositional (relative) abundances:
  

```r
p <- plot_composition(pseq2,
                      sample.sort = "nationality",
                      x.label = "nationality",
                      transform = "compositional") +
     guides(fill = guide_legend(ncol = 1))
```

```
## Error in plot_composition(pseq2, sample.sort = "nationality", x.label = "nationality", : object 'pseq2' not found
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```

```r
# Or try another theme
# from https://github.com/hrbrmstr/hrbrthemes
library(hrbrthemes)
library(gcookbook)
library(tidyverse)
```

```
## Loading tidyverse: tibble
## Loading tidyverse: tidyr
## Loading tidyverse: readr
## Loading tidyverse: purrr
```

```
## Warning: package 'tibble' was built under R version 3.3.3
```

```
## Warning: package 'readr' was built under R version 3.3.3
```

```
## Conflicts with tidy packages ----------------------------------------------
```

```
## filter(): dplyr, stats
## lag():    dplyr, stats
```

```r
p2 <- p + scale_y_percent() +
          labs(x="Samples", y="Relative abundance (%)",
                                   title="Relative abundance data",
                                   subtitle="Subtitle",
                                   caption="Caption here 'g'") + 
  theme_ipsum(grid="Y")
```

```
## Error in eval(expr, envir, enclos): object 'p' not found
```

```r
print(p2)  
```

```
## Error in print(p2): object 'p2' not found
```


Averaged by group:
  

```r
p <- plot_composition(pseq2,
                      average_by = "nationality", transform = "compositional")
```

```
## Error in plot_composition(pseq2, average_by = "nationality", transform = "compositional"): object 'pseq2' not found
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```



### Composition heatmaps


Plain heatmap on absolute abundances


```r
theme_set(theme_bw(30))
p <- plot_composition(pseq2, plot.type = "heatmap", mar = c(6, 13, 1, 1))
```


Heatmap with clr-transformed abundances:


```r
p <- plot_composition(pseq2, plot.type = "heatmap", transform = "clr", mar = c(6, 13, 1, 1))
```

```
## Error in plot_composition(pseq2, plot.type = "heatmap", transform = "clr", : object 'pseq2' not found
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```


Same with relative abundance, samples and OTUs sorted with the neatmap method:
  

```r
p <- plot_composition(pseq2, plot.type = "heatmap", transform = "compositional"
                      sample.sort = "neatmap", otu.sort = "neatmap",
                      mar = c(6, 13, 1, 1))
```


Same with Z-transformed, samples and OTUs sorted manually based on compositional data (Z-transformed data has negative values and the sorting method is not yet implemented for that):
  

```r
pseq3 <- transform(pseq2, "compositional")
sample.sort <- neatsort(pseq3, method = "NMDS", distance = "bray", target = "sites") 
otu.sort <- neatsort(pseq3, method = "NMDS", distance = "bray", target = "species")
p <- plot_composition(pseq2, plot.type = "heatmap", transform = "Z-OTU",
                      sample.sort = sample.sort, otu.sort = otu.sort,
                      mar = c(6, 13, 1, 1))
```

### Plot taxa prevalence

Here we use the Dynamics IBD data set from [Halfvarson J., et al. Nature Microbiology, 2017](http://www.nature.com/articles/nmicrobiol20174) as downloaded from [Qiita ID 1629](https://qiita.ucsd.edu/study/description/1629). 


```r
data(DynamicsIBD)
p0 <- DynamicsIBD
# Improve the taxonomic information
p0.f <- format_phyloseq(p0)
# For the available taxonomic levels, see tax_table(p0.f)
p <- plot_taxa_prevalence(p0.f, 'Phylum')
print(p)
```

![plot of chunk unnamed-chunk-1](figure/unnamed-chunk-1-1.png)

