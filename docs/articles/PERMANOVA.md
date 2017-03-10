---
title: "Comparisons"
author: "Leo Lahti"
date: "2017-03-05"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## PERMANOVA for community-level multivariate comparisons

PERMANOVA quantifies multivariate community-level differences between
groups.


Load example data:


```r
# Load libraries
library(microbiome)
library(ggplot2)
library(dplyr)

# Probiotics intervention example data 
data(peerj32) # Source: https://peerj.com/articles/32/
pseq <- peerj32$phyloseq # Rename the example data

# Pick relative abundances (compositional) and sample metadata 
pseq.rel <- transform(pseq, "compositional")
otu <- abundances(pseq.rel)
meta <- meta(pseq.rel)
```


### Visualize microbiome variation

Visualize the population density and highlight sample groups (probiotic treatment LGG  vs Placebo):


```r
p <- plot_landscape(pseq.rel, method = "NMDS", distance = "bray", col = "group", size = 3)
print(p)
```

<img src="figure/comparisons_permanova_visu-1.png" title="plot of chunk comparisons_permanova_visu" alt="plot of chunk comparisons_permanova_visu" width="300px" />


### PERMANOVA significance test for group-level differences

Now let us evaluate whether the group (probiotics vs. placebo) has a
significant effect on overall gut microbiota composition. Perform PERMANOVA: 


```r
# samples x species as input
library(vegan)
permanova <- adonis(t(otu) ~ group,
               data = meta, permutations=99, method = "bray")

# P-value
print(as.data.frame(permanova$aov.tab)["group", "Pr(>F)"])
```

```
## [1] 0.38
```


### Checking the homogeneity condition

Check that variance homogeneity assumptions hold (to ensure the reliability of the results):



```r
# Note the assumption of similar multivariate spread among the groups
# ie. analogous to variance homogeneity
# Here the groups have signif. different spreads and
# permanova result may be potentially explained by that.
dist <- vegdist(t(otu))
anova(betadisper(dist, meta$group))
```

```
## Analysis of Variance Table
## 
## Response: Distances
##           Df   Sum Sq   Mean Sq F value Pr(>F)
## Groups     1 0.000016 0.0000156  0.0042 0.9487
## Residuals 42 0.156962 0.0037372
```

### Investigate the top factors

Show coefficients for the top taxa separating the groups


```r
coef <- coefficients(permanova)["group1",]
top.coef <- coef[rev(order(abs(coef)))[1:20]]
par(mar = c(3, 14, 2, 1))
barplot(sort(top.coef), horiz = T, las = 1, main = "Top taxa")
```

<img src="figure/permanova_top-1.png" title="plot of chunk permanova_top" alt="plot of chunk permanova_top" width="300px" />

