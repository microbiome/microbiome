---
title: "microbiome vignette"
author: "Leo Lahti, Tineka Blake and Jarkko Salojarvi"
date: "2017-03-05"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteIndexEntry{microbiome tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}
-->




Introduction to the microbiome R package
===========

Tools for common analysis tasks in microbiome profiling studies;
illustrated with multiple example data sets from published studies;
extending the [phyloseq](http://joey711.github.io/phyloseq/import-data) class.


### Getting started

* [Installation](Template.md) 
* [Example data](Data.md)
* [Data manipulation](Preprocessing.md)


### Microbiome analysis

* [Alpha diversity](Diversity.md)
* [Beta diversity / Community heterogeneity](Betadiversity.md)
* [Community composition](Composition.md)
* [Core microbiota](Core.md)
* [Landscapes](Landscaping.md) (population density analysis)
* [Stability and tipping elements](Stability.md)


### Visualization and related tools

* [Heatmaps](Heatmap.md)
* [Networks](Networks.md)
* [Ordination](Ordination.md) (PCA, PCoA, NMDS, RDA etc.)
* [Regression](Regression.md)


### Statistical analysis

* [Bimodality](Bimodality.md)
* [Potential analysis](Potential.md)
* [Community comparisons](Comparisons.md) ([limma](limma.md), [PERMANOVA](PERMANOVA.md), [mixed models](Mixedmodels.md), [negative binomial](Negativebinomial.md)  etc)
* [Experimental tools](Experimental.md)




### Licensing and Citations

This work can be freely used, modified and distributed under the
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD\_licenses).

Kindly cite this work as follows:


```r
citation('microbiome')
```

```
## 
## To cite package 'microbiome' in publications use:
## 
##   Leo Lahti (2017). microbiome: Tools for microbiome analysis. R
##   package version 0.99.92. http://microbiome.github.com
## 
## A BibTeX entry for LaTeX users is
## 
##   @Manual{,
##     title = {microbiome: Tools for microbiome analysis},
##     author = {Leo Lahti},
##     year = {2017},
##     note = {R package version 0.99.92},
##     url = {http://microbiome.github.com},
##   }
```


### Dependencies

The package utilizes tools from a number of other R extensions,
including ade4 (Dray and Dufour, 2007; Chessel, Dufour, and Thioulouse, 2004; Dray, Dufour, and Chessel, 2007), dplyr (Wickham and Francois, 2016), ggplot2 (Wickham, 2009), MASS (Venables and Ripley, 2002), moments (Komsta and Novomestky, 2015), phyloseq (McMurdie and Holmes, 2013), RColorBrewer (Neuwirth, 2014), scales (Wickham, 2016), stats (R Core Team, 2016), tidyr (Wickham, 2017), vegan (Oksanen, Blanchet, Friendly, Kindt, Legendre, McGlinn, Minchin, O'Hara, Simpson, Solymos, Stevens, Szoecs, and Wagner, 2017).


### References



[1] D. Chessel, A. Dufour and J. Thioulouse. "The ade4 package-I-
One-table methods". In: _R News_ 4 (2004), pp. 5-10.

[2] S. Dray and A. Dufour. "The ade4 package: implementing the
duality diagram for ecologists". In: _Journal of Statistical
Software_ 22.4 (2007), pp. 1-20.

[3] S. Dray, A. Dufour and D. Chessel. "The ade4 package-II:
Two-table and K-table methods." In: _R News_ 7.2 (2007), pp.
47-52.

[4] L. Komsta and F. Novomestky. _moments: Moments, cumulants,
skewness, kurtosis and related tests_. R package version 0.14.
2015. <URL: https://CRAN.R-project.org/package=moments>.

[5] P. J. McMurdie and S. Holmes. "phyloseq: An R package for
reproducible interactive analysis and graphics of microbiome
census data". In: _PLoS ONE_ 8.4 (2013), p. e61217. <URL:
http://dx.plos.org/10.1371/journal.pone.0061217>.

[6] E. Neuwirth. _RColorBrewer: ColorBrewer Palettes_. R package
version 1.1-2. 2014. <URL:
https://CRAN.R-project.org/package=RColorBrewer>.

[7] J. Oksanen, F. G. Blanchet, M. Friendly, et al. _vegan:
Community Ecology Package_. R package version 2.4-2. 2017. <URL:
https://CRAN.R-project.org/package=vegan>.

[8] R Core Team. _R: A Language and Environment for Statistical
Computing_. R Foundation for Statistical Computing. Vienna,
Austria, 2016. <URL: https://www.R-project.org/>.

[9] W. N. Venables and B. D. Ripley. _Modern Applied Statistics
with S_. Fourth. ISBN 0-387-95457-0. New York: Springer, 2002.
<URL: http://www.stats.ox.ac.uk/pub/MASS4>.

[10] H. Wickham. _ggplot2: Elegant Graphics for Data Analysis_.
Springer-Verlag New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://ggplot2.org>.

[11] H. Wickham. _scales: Scale Functions for Visualization_. R
package version 0.4.1. 2016. <URL:
https://CRAN.R-project.org/package=scales>.

[12] H. Wickham. _tidyr: Easily Tidy Data with 'spread()' and
'gather()' Functions_. R package version 0.6.1. 2017. <URL:
https://CRAN.R-project.org/package=tidyr>.

[13] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.5.0. 2016. <URL:
https://CRAN.R-project.org/package=dplyr>.

