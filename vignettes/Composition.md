---
title: "Composition"
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
  %\VignetteIndexEntry{microbiome tutorial - composition}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Microbiota composition


### Barplots for composition

Also see [phyloseq barplot examples](http://joey711.github.io/phyloseq/plot_bar-examples.html) and [HITChip Barplots](Barplots.md)


Read example data from a [diet swap study](http://dx.doi.org/10.1038/ncomms7342):


```r
# Example data
library(microbiome)
library(ggplot2)
data(dietswap)
theme_set(theme_bw(22))
pseq <- dietswap
```

Show OTU absolute abundance in each sample. Plot absolute taxon
abundances (Note: on HITChip data the Phylum level is only
approximate):


```r
# Pick sample subset
library(phyloseq)
pseq2 <- subset_samples(pseq, group == "DI" & nationality == "AFR")
p <- plot_composition(pseq2, taxonomic.level = "Phylum") +
       theme(legend.position = "bottom") +
       guides(fill = guide_legend(nrow = 10, byrow = TRUE))
```

Arrange by sample variable and use custom X axis labels. Only consider the most abundant taxa. Africans have more Prevotella as expected. Absolute counts:


```r
# Pick the top OTUs only
top <- names(sort(taxa_sums(pseq), TRUE)[1:10])
pseq3 <- prune_taxa(top, pseq)
pseq3 <- subset_samples(pseq3, group == "DI" & timepoint.within.group == 1)

p <- plot_composition(pseq3, sample.sort = "nationality", x.label = "nationality") +
  guides(fill = guide_legend(ncol = 3)) +
  theme(legend.position = "bottom") +
  guides(fill = guide_legend(nrow = 5, byrow = TRUE))

print(p)
```

![plot of chunk composition-example4](figure/composition-example4-1.png)


Same with relative abundances:


```r
p <- plot_composition(pseq3, sample.sort = "nationality", x.label = "nationality", transformation = "relative.abundance") +
       guides(fill = guide_legend(ncol = 1))
print(p)
```

![plot of chunk composition-example4b](figure/composition-example4b-1.png)


Averaged by group:


```r
p <- plot_composition(pseq3, average_by = "nationality", transformation = "relative.abundance")
print(p)
```

![plot of chunk composition-example4c](figure/composition-example4c-1.png)



### Heatmaps for composition


Plain heatmap


```r
theme_set(theme_bw(30))
p <- plot_composition(pseq3, plot.type = "heatmap", mar = c(6, 13, 1, 1))
```

![plot of chunk composition-example5](figure/composition-example5-1.png)


Heatmap with Z-transformed OTUs


```r
p <- plot_composition(pseq3, plot.type = "heatmap", transformation = "Z-OTU", mar = c(6, 13, 1, 1))
```

![plot of chunk composition-example6](figure/composition-example6-1.png)


Same, but samples and OTUs sorted with the neatmap method:


```r
p <- plot_composition(pseq3, plot.type = "heatmap", transformation = "Z-OTU",
       			       sample.sort = "neatmap", otu.sort = "neatmap",
			       mar = c(6, 13, 1, 1))
```

```
## Warning in metaMDS(veganifyOTU(physeq), distance, ...): 'comm' has negative
## data: 'autotransform', 'noshare' and 'wascores' set to FALSE
```

```
## Warning in distfun(comm, method = distance, ...): results may be
## meaningless because data have negative entries in method "bray"
```

```
## Warning in metaMDS(veganifyOTU(physeq), distance, ...): 'comm' has negative
## data: 'autotransform', 'noshare' and 'wascores' set to FALSE
```

```
## Warning in distfun(comm, method = distance, ...): results may be
## meaningless because data have negative entries in method "bray"
```

```
## Warning in order_neatmap(x, method = "NMDS", distance = "bray", target = "species", : Order failed with species. 
##     			 Using default ordering.
```

![plot of chunk composition-example7](figure/composition-example7-1.png)


Same, but samples and OTUs sorted manually


```r
sample.sort <- order_neatmap(pseq3, method = "NMDS", distance = "bray", target = "sites", first = NULL) 
otu.sort <- order_neatmap(pseq3, method = "NMDS", distance = "bray", target = "species", first = NULL)

p <- plot_composition(pseq3, plot.type = "heatmap", transformation = "Z-OTU",
       			       sample.sort = sample.sort, otu.sort = otu.sort,
			       mar = c(6, 13, 1, 1))
```

![plot of chunk composition-example8](figure/composition-example8-1.png)



