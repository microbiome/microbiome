---
title: "Diversity"
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
  %\VignetteIndexEntry{microbiome tutorial - diversity}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Global Ecosystem State Variables 

Load example data:


```r
library(microbiome)
data(atlas1006)
pseq <- atlas1006
```


### Community richness and alpha diversity estimation

This function returns a table with a selection of global ecosystem indicators. Additional indices are also available; see the function help. See a separate page on [Beta diversity](Betadiversity.html).


```r
indicators <- global(pseq, measures = c("richness", "dominance"))
head(kable(indicators))
```

```
## [1] "|            | richness| dominance|"
## [2] "|:-----------|--------:|---------:|"
## [3] "|Sample-1    |      130| 0.1758679|"
## [4] "|Sample-2    |      130| 0.1716273|"
## [5] "|Sample-3    |      130| 0.2793253|"
## [6] "|Sample-4    |      130| 0.1957585|"
```

The supported divesity measures include those supported in the phyloseq::estimate_richness. Further measures are also provided (see function help), or can be calculated separately as described below.


### Dominance and rarity

The dominance and rarity indices refer to the abundance of the most and least abundant species, respectively. For the most abundant species, the sum of relative abundances up to the given rank (top-n) is calculated. For the least abundant taxa, the sum of relative abundances is returned for those taxa whose abundance falls below the indicated detection threshold.


```r
# Absolute abundances for the single most abundant taxa in each sample
ta <- dominance(pseq, rank = 1)

# Relative abundances
la <- rarity(pseq, detection = 0.2/100)
```


### Coverage

The coverage index gives the number of groups needed to have a given proportion of the ecosystem occupied (by default 0.5 ie 50%).


```r
do <- coverage(pseq, threshold = 0.5)
```


### Rarity and core abundance

The core_abundance function refers to the relative proportion of the core species. Rare abundance provides the complement.


```r
ra <- rare_abundance(pseq, detection = .1/100, prevalence = 50/100)
co <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)
```



### Gini index

Gini index is a common measure for inequality in economical income, but can also be used as a community diversity measure.


```r
gi <- inequality(pseq)
```




### Visualization

Show indicators:


```r
library(ggplot2)
theme_set(theme_bw(20)) # Set bw color scheme
p <- ggplot(indicators, aes(x = richness)) + geom_histogram()
print(p)
```

![plot of chunk div-example2](figure/div-example2-1.png)

### Group-wise comparison

Visualize ecosystem state indicators w.r.t. discrete variable (or check more generic [group-wise comparison tools](Comparisons.html)). You can also indicate subjects by lines (assuming the subject field is available in the metadata):


```r
p <- plot_diversity(pseq, "bmi_group", measures = c("chao1", "shannon"), indicate.subjects = TRUE)
print(p)
```

![plot of chunk div-example2bb](figure/div-example2bb-1.png)

Indicators vs. continuous variable:


```r
library(dplyr)
pseq <- atlas1006

# Add shannon diversity in the sample metadata
sample_data(pseq)$diversity <- global(pseq, measures = "shannon")$shannon

# Visualize
p <- plot_regression(diversity ~ age, meta(pseq))
print(p)
```

![plot of chunk indicators-example13](figure/indicators-example13-1.png)


