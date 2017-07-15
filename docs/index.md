<!--
  %\VignetteIndexEntry{microbiome tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}
-->
[![Build
Status](https://api.travis-ci.org/microbiome/microbiome.png)](https://travis-ci.org/microbiome/microbiome)

Tools for microbiome analysis; with multiple example data sets from
published studies; extending the
[phyloseq](http://joey711.github.io/phyloseq/import-data) class.

Installation and use
====================

Getting started
---------------

-   [Overview (vignette)](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)
-   [Installation](Template.html)
-   [Example data](Data.html)
-   [Data manipulation](Preprocessing.html)

Microbiome analysis
-------------------

-   [Alpha diversity](Diversity.html)
-   [Beta diversity / Community heterogeneity](Betadiversity.html)
-   [Community composition](Composition.html)
-   [Core microbiota](Core.html)
-   [Landscapes](Landscaping.html) (population density analysis)
-   [Stability and tipping elements](Stability.html)

Visualization and related tools
-------------------------------

-   [Heatmaps](Heatmap.html)
-   [Networks](Networks.html)
-   [Ordination](Ordination.html) (PCA, PCoA, NMDS, RDA etc.)
-   [Regression](Regression.html)

Statistical analysis
--------------------

-   [Bimodality](Bimodality.html)
-   [Community comparisons](Comparisons.html) ([limma](limma.html),
    [PERMANOVA](PERMANOVA.html), [mixed models](Mixedmodels.html),
    [negative binomial](Negativebinomial.html) etc)
-   [Experimental tools](Experimental.html)

Development
===========

New examples, tutorial pages, and other contributions are
[welcome](Contributing.html). The material can be freely used, modified
and distributed under the [Two-clause FreeBSD
license](http://en.wikipedia.org/wiki/BSD_licenses). For source code,
see the [Github page](https://github.com/microbiome/microbiome/).

Acknowledgements
================

**Kindly cite this work** as follows: "Leo Lahti, Sudarshan Shetty [et
al.](https://github.com/microbiome/microbiome/graphs/contributors)
(2017). Tools for microbiome analysis in R. Version 0.99.49. URL:
<http://microbiome.github.com/microbiome>. Check also the relevant
references listed in the manual page of each function.

The package utilizes tools from a number of other R extensions,
including dplyr (Wickham, Francois, Henry, and Müller, 2017), ggplot2
(Wickham, 2009), phyloseq (McMurdie and Holmes, 2013), tidyr (Wickham,
2017), vegan (Oksanen, Blanchet, Friendly, Kindt, Legendre, McGlinn,
Minchin, O'Hara, Simpson, Solymos, Stevens, Szoecs, and Wagner, 2017).
