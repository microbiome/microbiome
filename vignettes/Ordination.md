---
title: "Ordination"
author: "Leo Lahti"
date: "2016-11-03"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - ordination}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Ordination examples

Some examples with HITChip data. See also [phyloseq ordination tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html).

Load example data:


```r
library(microbiome)
library(phyloseq)
library(ggplot2)
data(atlas1006)
pseq <- atlas1006

# Convert signal to relative abundances
pseq.rel <- transform_phyloseq(pseq, "relative.abundance")

# Pick OTUs that are present with >1 percent relative abundance 
# in >10 percent of the samples
pseq2 <- filter_taxa(pseq.rel, function(x) sum(x > 1) > (0.1*nsamples(pseq.rel)), TRUE)
```


### Sample ordination

Project the samples with the given method and distance. See also plot_ordination from the phyloseq package.



```r
set.seed(423542)
pseq.ord <- ordinate(pseq2, "NMDS", "bray")
# Just pick the projected data (first two columns + metadata)
proj <- plot_ordination(pseq2, pseq.ord, justDF = T)
```

Then visualize the projected data:


```r
# Highlighting gender
p <- densityplot(proj[, 1:2], col = proj$gender, legend = T)
```

```
## Error in UseMethod("densityplot"): no applicable method for 'densityplot' applied to an object of class "data.frame"
```

```r
print(p)

# Projection with sample names:
ax1 <- names(proj)[[1]]
ax2 <- names(proj)[[2]]
p <- ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) +
       geom_text(size = 2)
print(p)
```

<img src="figure/ordination4-1.png" title="plot of chunk ordination4" alt="plot of chunk ordination4" width="400px" /><img src="figure/ordination4-2.png" title="plot of chunk ordination4" alt="plot of chunk ordination4" width="400px" />


Ordinate the taxa in NMDS plot with Bray-Curtis distances


```r
p <- plot_ordination(pseq2, pseq.ord, type = "taxa", color = "Phylum", title = "Taxa ordination")
print(p)
```

![plot of chunk ordination-pca-ordination21](figure/ordination-pca-ordination21-1.png)

Grouping by Phyla


```r
p + facet_wrap(~Phylum, 5)
```

![plot of chunk ordination-pca-ordination22](figure/ordination-pca-ordination22-1.png)


### Multidimensional scaling (MDS / PCoA)


```r
plot_ordination(pseq, ordinate(pseq, "MDS"), color = "DNA_extraction_method") + geom_point(size = 5)
```

![plot of chunk ordination-ordinate23](figure/ordination-ordinate23-1.png)

### RDA

See a separate page on [RDA](RDA.md).


### Canonical correspondence analysis (CCA)



```r
# With samples
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "samples", color = "gender")
p <- p + geom_point(size = 4)
print(p)

# With taxa:
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "taxa", color = "Phylum")
p <- p + geom_point(size = 4)
print(p)
```

<img src="figure/ordination-ordinate24a-1.png" title="plot of chunk ordination-ordinate24a" alt="plot of chunk ordination-ordinate24a" width="400px" /><img src="figure/ordination-ordinate24a-2.png" title="plot of chunk ordination-ordinate24a" alt="plot of chunk ordination-ordinate24a" width="400px" />


### Split plot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "split", shape = "gender", 
    color = "Phylum", label = "gender")
```

![plot of chunk ordination-ordinate25](figure/ordination-ordinate25-1.png)


### Ordination biplot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "biplot", color = "Phylum")
```

![plot of chunk ordination-ordinate26](figure/ordination-ordinate26-1.png)






