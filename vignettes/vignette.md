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
 * scales (Wickham, 2016b)
 * stats (R Core Team, 2016)
 * tidyr (Wickham, 2016a)
 * vegan (Oksanen, Blanchet, Kindt, et al., 2016)


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
##  [1] splines   stats4    tcltk     parallel  grid      stats     graphics 
##  [8] grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] googleVis_0.5.10      scales_0.4.0          SpiecEasi_0.1        
##  [4] VGAM_1.0-1            huge_1.2.7            igraph_1.0.1         
##  [7] Matrix_1.2-5          earlywarnings_1.1.22  tseries_0.10-34      
## [10] tgp_2.4-14            moments_0.14          ggrepel_0.5.1        
## [13] gridExtra_2.2.1       FD_1.0-12             geometry_0.3-6       
## [16] magic_1.5-6           abind_1.4-3           ape_3.4              
## [19] ade4_1.7-4            vegan_2.3-5           lattice_0.20-33      
## [22] permute_0.9-0         HITChipDB_0.6.30      RPA_1.27.41          
## [25] affy_1.48.0           Biobase_2.30.0        BiocGenerics_0.16.1  
## [28] RMySQL_0.10.8         preprocessCore_1.32.0 knitcitations_1.0.7  
## [31] knitr_1.12.3          intergraph_2.0-2      sna_2.3-2            
## [34] network_1.13.0        ggnet_0.1.0           GGally_1.0.1         
## [37] devtools_1.11.0       limma_3.26.9          sorvi_0.7.45         
## [40] tibble_1.0            ggplot2_2.1.0         tidyr_0.4.1          
## [43] dplyr_0.4.3           MASS_7.3-45           netresponse_1.20.15  
## [46] reshape2_1.4.1        mclust_5.2            minet_3.28.0         
## [49] Rgraphviz_2.14.0      graph_1.48.0          microbiome_0.99.81   
## [52] RSQLite_1.0.0         DBI_0.3.1             phyloseq_1.14.0      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      dynamicTreeCut_1.63-1 som_0.3-5            
##  [4] qvalue_2.2.2          XVector_0.10.0        affyio_1.40.0        
##  [7] AnnotationDbi_1.32.3  mvtnorm_1.0-5         lubridate_1.5.6      
## [10] RefManageR_0.10.13    codetools_0.2-14      doParallel_1.0.10    
## [13] impute_1.44.0         spam_1.3-0            Formula_1.2-1        
## [16] Cairo_1.5-9           WGCNA_1.51            cluster_2.0.4        
## [19] GO.db_3.2.2           Kendall_2.2           httr_1.1.0           
## [22] lazyeval_0.1.10       assertthat_0.1        formatR_1.3          
## [25] acepack_1.3-3.3       tools_3.2.5           gtable_0.2.0         
## [28] maps_3.1.0            Rcpp_0.12.4           Biostrings_2.38.4    
## [31] RJSONIO_1.3-0         multtest_2.26.0       biom_0.3.12          
## [34] nlme_3.1-127          iterators_1.0.8       lmtest_0.9-34        
## [37] fastcluster_1.1.20    stringr_1.0.0         XML_3.98-1.4         
## [40] zlibbioc_1.16.0       zoo_1.7-12            BiocInstaller_1.20.1 
## [43] RColorBrewer_1.1-2    fields_8.3-6          curl_0.9.7           
## [46] memoise_1.0.0         rpart_4.1-10          reshape_0.8.5        
## [49] latticeExtra_0.6-28   stringi_1.0-1         maptree_1.4-7        
## [52] highr_0.5.1           S4Vectors_0.8.11      foreach_1.4.3        
## [55] nortest_1.0-4         boot_1.3-18           bibtex_0.4.0         
## [58] chron_2.3-47          matrixStats_0.50.2    bitops_1.0-6         
## [61] dmt_0.8.20            evaluate_0.8.3        labeling_0.3         
## [64] plyr_1.8.3            magrittr_1.5          R6_2.1.2             
## [67] IRanges_2.4.8         Hmisc_3.17-3          foreign_0.8-66       
## [70] withr_1.0.1           mgcv_1.8-12           survival_2.39-2      
## [73] RCurl_1.95-4.8        nnet_7.3-12           KernSmooth_2.23-15   
## [76] data.table_1.9.6      git2r_0.14.0          digest_0.6.9         
## [79] munsell_0.4.3         quadprog_1.5-5
```




