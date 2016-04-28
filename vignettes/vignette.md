---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2016-04-28"
bibliography: 
- bibliography.bib
- references.bib
output: 
  md_document:
    variant: markdown_github
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial}
  %\usepackage[utf8]{inputenc}
-->




microbiome R package
===========

Tools for microbiome analysis in R (R Core Team, 2013).

The microbiome R package provides analysis tools, diagnostic plots,
and other utilities for microbiome analyses, as well as multiple
example data sets.

We use the [phyloseq](http://joey711.github.io/phyloseq/import-data)
class, a standard representation format in R for taxonomic
profiling. This package provides extensions and convenient wrappers
for many standard tasks encountered in microbiome studies. 


### Getting started

* [Installation](Template.md) 
* [Example data](Data.md)
* [Data manipulation](Preprocessing.md)


### Visualization and related tools

* [Community composition](Composition.md)
* [Heatmaps](Heatmap.md)
* [Networks](Networks.md)
* [Ordination (PCA, PCoA, NMDS etc.)](Ordination.md)
* [Density](Density.md)
* [Interactive](Interactive.md)


### Microbiota composition

* [Bistability analysis](Stability.md)
* [Variability](Variability.md) (Intra- and inter-individual 'stability')
* [Core microbiota](Core.md)
* [Diversity](Diversity.md)


### Analysis

* [Pairwise comparisons](Comparisons.md)
* [Stability, bimodality, tipping elements](Stability.md)
* [Clustering](Clustering.md)
* [Linear models](limma.md)
* [Phylogenetic microarrays](HITChip.Rmd)


### Licensing and Citations

This work can be freely used, modified and distributed under the
[Two-clause FreeBSD
license](http://en.wikipedia.org/wiki/BSD\_licenses). Kindly cite as
'Leo Lahti and Jarkko Salojarvi (2014-2016). microbiome R
package. URL: http://microbiome.github.com'.


### Dependencies

The package utilizes tools from a number of other R extensions,
including:

 * ade4 (Dray and Dufour, 2007; Chessel, Dufour, and Thioulouse, 2004; Dray, Dufour, and Chessel, 2007)
 * dplyr (Wickham and Francois, 2015)  
 * ggplot2 (Wickham, 2009)
 * MASS (Venables and Ripley, 2002)
 * moments (Komsta and Novomestky, 2015)
 * netresponse (Lahti, Knuuttila, and Kaski, 2010) 
 * phyloseq (McMurdie and Holmes, 2013)
 * RColorBrewer (Neuwirth, 2014)
 * scales (Wickham, 2016)
 * stats (R Core Team, 2016)
 * tidyr (Wickham, 2016)
 * vegan (Oksanen, Blanchet, Kindt, Legendre, Minchin, O'Hara, Simpson, Solymos, Stevens, and Wagner, 2016)


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

[5] L. Lahti, J. E. Knuuttila and S. Kaski. "Global modeling of
transcriptional responses in interaction networks". In:
_Bioinformatics_ 26 (21 2010), pp. 2713-20.

[6] P. J. McMurdie and S. Holmes. "phyloseq: An R package for
reproducible interactive analysis and graphics of microbiome
census data". In: _PLoS ONE_ 8.4 (2013), p. e61217. <URL:
http://dx.plos.org/10.1371/journal.pone.0061217>.

[7] E. Neuwirth. _RColorBrewer: ColorBrewer Palettes_. R package
version 1.1-2. 2014. <URL:
https://CRAN.R-project.org/package=RColorBrewer>.

[8] J. Oksanen, F. G. Blanchet, R. Kindt, et al. _vegan: Community
Ecology Package_. R package version 2.3-5. 2016. <URL:
https://CRAN.R-project.org/package=vegan>.

[9] R Core Team. _R: A language and environment for statistical
computing_. Vienna, Austria: R Foundation for Statistical
Computing, 2013. ISBN: ISBN 3-900051-07-0. <URL:
http://www.R-project.org/>.

[10] R Core Team. _R: A Language and Environment for Statistical
Computing_. R Foundation for Statistical Computing. Vienna,
Austria, 2016. <URL: https://www.R-project.org/>.

[11] W. N. Venables and B. D. Ripley. _Modern Applied Statistics
with S_. Fourth. ISBN 0-387-95457-0. New York: Springer, 2002.
<URL: http://www.stats.ox.ac.uk/pub/MASS4>.

[12] H. Wickham. _ggplot2: Elegant Graphics for Data Analysis_.
Springer-Verlag New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://ggplot2.org>.

[13] H. Wickham. _scales: Scale Functions for Visualization_. R
package version 0.4.0. 2016. <URL:
https://CRAN.R-project.org/package=scales>.

[14] H. Wickham. _tidyr: Easily Tidy Data with `spread()` and
`gather()` Functions_. R package version 0.4.1. 2016. <URL:
https://CRAN.R-project.org/package=tidyr>.

[15] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.4.3. 2015. <URL:
https://CRAN.R-project.org/package=dplyr>.

### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.5 (2016-04-14)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_BE.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=de_BE.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=de_BE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_BE.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] tcltk     parallel  grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] FD_1.0-12             geometry_0.3-6        magic_1.5-6          
##  [4] abind_1.4-3           ape_3.4               ade4_1.7-4           
##  [7] vegan_2.3-5           lattice_0.20-33       permute_0.9-0        
## [10] HITChipDB_0.6.30      RPA_1.27.41           affy_1.48.0          
## [13] Biobase_2.30.0        BiocGenerics_0.16.1   RMySQL_0.10.8        
## [16] preprocessCore_1.32.0 knitcitations_1.0.7   knitr_1.12.3         
## [19] intergraph_2.0-2      sna_2.3-2             network_1.13.0       
## [22] ggnet_0.1.0           GGally_1.0.1          devtools_1.11.0      
## [25] limma_3.26.9          sorvi_0.7.45          tibble_1.0           
## [28] ggplot2_2.1.0         tidyr_0.4.1           dplyr_0.4.3          
## [31] MASS_7.3-45           netresponse_1.20.15   reshape2_1.4.1       
## [34] mclust_5.2            minet_3.28.0          Rgraphviz_2.14.0     
## [37] graph_1.48.0          microbiome_0.99.81    RSQLite_1.0.0        
## [40] DBI_0.3.1             phyloseq_1.14.0      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      dynamicTreeCut_1.63-1 som_0.3-5            
##  [4] qvalue_2.2.2          XVector_0.10.0        affyio_1.40.0        
##  [7] AnnotationDbi_1.32.3  mvtnorm_1.0-5         lubridate_1.5.6      
## [10] RefManageR_0.10.13    codetools_0.2-14      splines_3.2.5        
## [13] doParallel_1.0.10     impute_1.44.0         tgp_2.4-14           
## [16] spam_1.3-0            Formula_1.2-1         WGCNA_1.51           
## [19] cluster_2.0.4         GO.db_3.2.2           Kendall_2.2          
## [22] httr_1.1.0            lazyeval_0.1.10       assertthat_0.1       
## [25] Matrix_1.2-5          formatR_1.3           acepack_1.3-3.3      
## [28] tools_3.2.5           igraph_1.0.1          gtable_0.2.0         
## [31] maps_3.1.0            Rcpp_0.12.4           Biostrings_2.38.4    
## [34] RJSONIO_1.3-0         multtest_2.26.0       biom_0.3.12          
## [37] nlme_3.1-127          iterators_1.0.8       lmtest_0.9-34        
## [40] fastcluster_1.1.20    stringr_1.0.0         XML_3.98-1.4         
## [43] zlibbioc_1.16.0       zoo_1.7-12            scales_0.4.0         
## [46] BiocInstaller_1.20.1  RColorBrewer_1.1-2    fields_8.3-6         
## [49] memoise_1.0.0         gridExtra_2.2.1       rpart_4.1-10         
## [52] reshape_0.8.5         latticeExtra_0.6-28   stringi_1.0-1        
## [55] maptree_1.4-7         highr_0.5.1           S4Vectors_0.8.11     
## [58] tseries_0.10-34       foreach_1.4.3         nortest_1.0-4        
## [61] boot_1.3-18           bibtex_0.4.0          chron_2.3-47         
## [64] moments_0.14          matrixStats_0.50.1    bitops_1.0-6         
## [67] dmt_0.8.20            evaluate_0.8.3        labeling_0.3         
## [70] plyr_1.8.3            magrittr_1.5          R6_2.1.2             
## [73] IRanges_2.4.8         earlywarnings_1.1.22  Hmisc_3.17-3         
## [76] foreign_0.8-66        withr_1.0.1           mgcv_1.8-12          
## [79] survival_2.39-2       RCurl_1.95-4.8        nnet_7.3-12          
## [82] KernSmooth_2.23-15    data.table_1.9.6      digest_0.6.9         
## [85] stats4_3.2.5          munsell_0.4.3         quadprog_1.5-5
```




