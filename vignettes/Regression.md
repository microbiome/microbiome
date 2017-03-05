---
title: "Regression"
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
  %\VignetteIndexEntry{microbiome tutorial - regression}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


## Regression plots

Load example data


```r
library(microbiome)
data(atlas1006)
pseq <- atlas1006
```

Regression curve with smoothed error bars based on
the [Visually-Weighted Regression](http://www.fight-entropy.com/2012/07/visually-weighted-regression.html) by Solomon M. Hsiang. The sorvi implementation extends [Felix Schonbrodt's original code](http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/). See also [potential analysis](Potential.md).


```r
p <- plot_regression(diversity ~ age, meta(pseq))
print(p)
```

![plot of chunk variability-regression](figure/variability-regression-1.png)

