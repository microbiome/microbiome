<!--
  %\VignetteIndexEntry{microbiome tutorial}
  %\VignetteEngine{knitr::rmarkdown}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}
-->
Introduction to the microbiome R package
========================================

[![Build
Status](https://api.travis-ci.org/microbiome/microbiome.png)](https://travis-ci.org/microbiome/microbiome)

Tools for common analysis tasks in microbiome profiling studies;
illustrated with multiple example data sets from published studies;
extending the [phyloseq](http://joey711.github.io/phyloseq/import-data)
class.

### Getting started

-   [Overview (vignette)](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)
-   [Installation](Template.html)
-   [Example data](Data.html)
-   [Data manipulation](Preprocessing.html)

### Microbiome analysis

-   [Alpha diversity](Diversity.html)
-   [Beta diversity / Community heterogeneity](Betadiversity.html)
-   [Community composition](Composition.html)
-   [Core microbiota](Core.html)
-   [Landscapes](Landscaping.html) (population density analysis)
-   [Stability and tipping elements](Stability.html)

### Visualization and related tools

-   [Heatmaps](Heatmap.html)
-   [Networks](Networks.html)
-   [Ordination](Ordination.html) (PCA, PCoA, NMDS, RDA etc.)
-   [Regression](Regression.html)

### Statistical analysis

-   [Bimodality](Bimodality.html)
-   [Community comparisons](Comparisons.html) ([limma](limma.html),
    [PERMANOVA](PERMANOVA.html), [mixed models](Mixedmodels.html),
    [negative binomial](Negativebinomial.html) etc)
-   [Experimental tools](Experimental.html)

### Development

New examples and tutorial pages from the user community are welcome:

-   [Contributing](Contributing.html)

### Licensing and Citations

This work can be freely used, modified and distributed under the
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD_licenses).

Kindly cite this work as follows:

    citation('microbiome')

    ## 
    ## Please cite the package as follows:
    ## 
    ##   Leo Lahti et al.  microbiome R package.  URL:
    ##   http://microbiome.github.io
    ## 
    ## A BibTeX entry for LaTeX users is
    ## 
    ##   @Misc{,
    ##     title = {microbiome R package},
    ##     author = {Leo Lahti and Sudarshan Shetty and Tineka Blake and Jarkko Salojarvi},
    ##     year = {2012-2017},
    ##   }

### Dependencies

The package utilizes tools from a number of other R extensions,
including ade4 (Dray and Dufour, 2007; Chessel, Dufour, and Thioulouse,
2004; Dray, Dufour, and Chessel, 2007), dplyr (Wickham, Francois, Henry,
and Müller, 2017), ggplot2 (Wickham, 2009), MASS (Venables and Ripley,
2002), moments (Komsta and Novomestky, 2015), phyloseq (McMurdie and
Holmes, 2013), RColorBrewer (Neuwirth, 2014), scales (Wickham, 2016),
stats (R Core Team, 2017), tidyr (Wickham, 2017), vegan (Oksanen,
Blanchet, Friendly, Kindt, Legendre, McGlinn, Minchin, O'Hara, Simpson,
Solymos, Stevens, Szoecs, and Wagner, 2017).

### References

\[1\] D. Chessel, A. Dufour and J. Thioulouse. "The ade4 package-I-
One-table methods". In: *R News* 4 (2004), pp. 5-10. \[1\] S. Dray and
A. Dufour. "The ade4 package: implementing the duality diagram for
ecologists". In: *Journal of Statistical Software* 22.4 (2007), pp.
1-20. \[1\] S. Dray, A. Dufour and D. Chessel. "The ade4 package-II:
Two-table and K-table methods." In: *R News* 7.2 (2007), pp. 47-52.
\[1\] L. Komsta and F. Novomestky. *moments: Moments, cumulants,
skewness, kurtosis and related tests*. R package version 0.14. 2015.
&lt;URL: <https://CRAN.R-project.org/package=moments>&gt;. \[1\] P. J.
McMurdie and S. Holmes. "phyloseq: An R package for reproducible
interactive analysis and graphics of microbiome census data". In: *PLoS
ONE* 8.4 (2013), p. e61217. &lt;URL:
<http://dx.plos.org/10.1371/journal.pone.0061217>&gt;. \[1\] E.
Neuwirth. *RColorBrewer: ColorBrewer Palettes*. R package version 1.1-2.
2014. &lt;URL: <https://CRAN.R-project.org/package=RColorBrewer>&gt;.
\[1\] J. Oksanen, F. G. Blanchet, M. Friendly, et al. *vegan: Community
Ecology Package*. R package version 2.4-3. 2017. &lt;URL:
<https://CRAN.R-project.org/package=vegan>&gt;. \[1\] R Core Team. *R: A
Language and Environment for Statistical Computing*. R Foundation for
Statistical Computing. Vienna, Austria, 2017. &lt;URL:
<https://www.R-project.org/>&gt;. \[1\] W. N. Venables and B. D. Ripley.
*Modern Applied Statistics with S*. Fourth. ISBN 0-387-95457-0. New
York: Springer, 2002. &lt;URL:
<http://www.stats.ox.ac.uk/pub/MASS4>&gt;. \[1\] H. Wickham. *ggplot2:
Elegant Graphics for Data Analysis*. Springer-Verlag New York, 2009.
ISBN: 978-0-387-98140-6. &lt;URL: <http://ggplot2.org>&gt;. \[1\] H.
Wickham. *scales: Scale Functions for Visualization*. R package version
0.4.1. 2016. &lt;URL: <https://CRAN.R-project.org/package=scales>&gt;.
\[1\] H. Wickham. *tidyr: Easily Tidy Data with 'spread()' and
'gather()' Functions*. R package version 0.6.3. 2017. &lt;URL:
<https://CRAN.R-project.org/package=tidyr>&gt;. \[1\] H. Wickham, R.
Francois, L. Henry, et al. *dplyr: A Grammar of Data Manipulation*. R
package version 0.7.0. 2017. &lt;URL:
<https://CRAN.R-project.org/package=dplyr>&gt;.
