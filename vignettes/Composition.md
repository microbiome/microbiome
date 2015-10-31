## Microbiota composition


### Barplot visualizations

Also see [phyloseq barplot examples](http://joey711.github.io/phyloseq/plot_bar-examples.html) and [HITChip Barplots](Barplots.md)

Loading example data:


```r
# Example data
library(microbiome)
```

```
## Loading required package: phyloseq
## Creating a generic function for 'nchar' from package 'base' in package 'S4Vectors'
## Loading required package: RPA
## Loading required package: parallel
## Loading required package: affy
## Loading required package: BiocGenerics
## 
## Attaching package: 'BiocGenerics'
## 
## The following objects are masked from 'package:parallel':
## 
##     clusterApply, clusterApplyLB, clusterCall, clusterEvalQ,
##     clusterExport, clusterMap, parApply, parCapply, parLapply,
##     parLapplyLB, parRapply, parSapply, parSapplyLB
## 
## The following object is masked from 'package:stats':
## 
##     xtabs
## 
## The following objects are masked from 'package:base':
## 
##     anyDuplicated, append, as.data.frame, as.vector, cbind,
##     colnames, do.call, duplicated, eval, evalq, Filter, Find, get,
##     intersect, is.unsorted, lapply, Map, mapply, match, mget,
##     order, paste, pmax, pmax.int, pmin, pmin.int, Position, rank,
##     rbind, Reduce, rep.int, rownames, sapply, setdiff, sort,
##     table, tapply, union, unique, unlist, unsplit
## 
## Loading required package: Biobase
## Welcome to Bioconductor
## 
##     Vignettes contain introductory material; view with
##     'browseVignettes()'. To cite Bioconductor, see
##     'citation("Biobase")', and for packages 'citation("pkgname")'.
## 
## 
## Attaching package: 'Biobase'
## 
## The following object is masked from 'package:phyloseq':
## 
##     sampleNames
## 
## 
## RPA Copyright (C) 2008-2014 Leo Lahti.
## This program comes with ABSOLUTELY NO WARRANTY.
## This is free software, and you are welcome to redistribute it under the FreeBSD open source license.
## 
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
```

```r
pseq0 <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```r
# Pick sample subset
pseq <- subset_samples(pseq0, group == "DI" & nationality == "AFR")
```

Show OTU absolute abundance in each sample. Plot absolute taxon
abundances (Note: on HITChip data the Phylum level is only
approximate):


```r
plot_composition(pseq, taxonomic.level = "Phylum")
```

![plot of chunk composition-example1b](figure/composition-example1b-1.png) 

Same with relative abundances:


```r
p <- plot_composition(pseq, taxonomic.level = "Phylum", relative.abundance = TRUE)
p <- p + guides(fill = guide_legend(nrow = 12, byrow = TRUE))
```

```
## Error in eval(expr, envir, enclos): could not find function "guides"
```

```r
p <- p + theme(legend.position = "bottom")
```

```
## Error in eval(expr, envir, enclos): could not find function "theme"
```

```r
print(p)
```

![plot of chunk composition-example3](figure/composition-example3-1.png) 


Arrange by sample variable and use custom X axis labels. Africans have more Prevotella as expected:


```r
# Subset taxa and samples
pseq <- subset_samples(pseq0, group == "DI" & timepoint.within.group == 1)
# Pick the top OTUs only
pseq <- prune_taxa(names(sort(taxa_sums(pseq), TRUE)[1:5]), pseq)
p <- plot_composition(pseq, relative.abundance = TRUE, sort.by = "nationality", x.label = "nationality")
p <- p + guides(fill = guide_legend(ncol = 1))
```

```
## Error in eval(expr, envir, enclos): could not find function "guides"
```

```r
p <- p + theme(legend.position = "bottom")
```

```
## Error in eval(expr, envir, enclos): could not find function "theme"
```

```r
print(p)
```

![plot of chunk composition-example4](figure/composition-example4-1.png) 

### Coloured Barplots

The following example visualizes samples, colored by Phylum
percentages (in this example data the Phylum is approximated by 16S
sequence similarity, not exactly Phylum):


```r
pseq <- subset_samples(pseq0, group == "DI")
p <- plot_bar(pseq, x = "timepoint.within.group", fill = "Phylum", facet_grid = ~nationality)
print(p)
```

![plot of chunk composition-example5](figure/composition-example5-1.png) 

