---
title: "RDA"
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
  %\VignetteIndexEntry{microbiome tutorial - rda}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## RDA analysis and visualization. 

Load the package and example data:


```r
library(microbiome)
# Data from https://peerj.com/articles/32/
data("peerj32")
pseq <- peerj32$phyloseq
```

### Standard RDA 

Standard RDA for microbiota profiles versus the given (here 'time')
variable from sample metadata:


```r
pseq.trans <- transform_phyloseq(pseq, "hell") # Hellinger transformation
rda.result <- rda_physeq(pseq.trans, "time", scale = TRUE)

# Proportion explained by the contraints
summary(rda.result)$constr.chi/summary(rda.result)$tot.chi
```

```
## [1] 0.01540884
```

### RDA visualization

Visualizing the standard RDA output:


```r
library(phyloseq)
meta <- sample_data(pseq.trans)
plot(rda.result, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = meta$time)
points(rda.result, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = meta$time)
library(vegan)
pl <- ordihull(rda.result, meta$time, scaling = 3, label = TRUE)
title("RDA")
```

![plot of chunk rda4](figure/rda4-1.png)

See also the RDA method in phyloseq::ordinate, which is calculated without the formula.


### RDA significance test


```r
library(vegan)
permutest(rda.result) 
```

```
## 
## Permutation test for rda 
## 
## Permutation: free
## Number of permutations: 99
##  
## Call: rda(formula = otu ~ annot, scale = scale, na.action =
## na.action)
## Permutation test for all constrained eigenvalues
## Pseudo-F:	 0.6572996 (with 1, 42 Degrees of Freedom)
## Significance:	 0.85
```

### Bagged RDA

Fitting bagged (bootstrap aggregated) RDA on a phyloseq object:


```r
res <- bagged_rda(pseq.trans, "group", sig.thresh=0.05, nboot=100)
```

Visualizing bagged RDA:


```r
plot_bagged_rda(res)
```

![plot of chunk rda6](figure/rda6-1.png)


### RDA with confounding variables 

For more complex RDA scenarios, use the vegan package:


```r
# Pick microbiota profiling data from the phyloseq object
otu <- taxa_abundances(pseq.trans)

# Sample annotations
meta <- sample_data(pseq.trans)

# RDA with confounders
rda.result2 <- rda(t(otu) ~ meta$time + Condition(meta$subject + meta$gender))
```



