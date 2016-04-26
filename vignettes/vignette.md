---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2016-04-26"
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

Tools for microbiome analysis in R ((R Core Team, 2013)).


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

The package utilizes tools from a number of other CRAN and
Bioconductor extensions, including:

 * ade4 (Dray and Dufour, 2007; Chessel, Dufour, and Thioulouse, 2004; Dray, Dufour, and Chessel, 2007)
 * df2json (Caballero, 2013)
 * dplyr (Wickham and Francois, 2015)  
 * fastcluster (Müllner, 2013)
 * ggplot2 (Wickham, 2009)
 * minet (Meyer, Lafitte, and Bontempi, 2008) 
 * mixOmics (Cao, Gonzalez, Rohart, Gautier, Monget, Coquery, Yao, and Liquet., 2015)
 * phyloseq (McMurdie and Holmes, 2013)
 * RPA (Lahti, Torrente, Elo, Brazma, and Rung, 2013; Lahti, Elo, Aittokallio, and Kaski, 2011) 
 * reshape (Wickham and Hadley, 2007) 
 * rjson (Couture-Beil, 2014)
 * RCurl (Temple Lang and team, 2016)
 * vegan (Oksanen, Blanchet, Kindt, Legendre, Minchin, O'Hara, Simpson, Solymos, Stevens, and Wagner, 2016)
 * WGCNA (Langfelder and Horvath, 2008; Langfelder and Horvath, 2012)


### References



[1] N. Caballero. _df2json: Convert a dataframe to JSON_. R
package version 0.0.2. 2013. <URL:
https://CRAN.R-project.org/package=df2json>.

[2] K. L. Cao, I. Gonzalez, S. D. w. k. c. F. Rohart, et al.
_mixOmics: Omics Data Integration Project_. R package version
5.2.0. 2015. <URL: https://CRAN.R-project.org/package=mixOmics>.

[3] D. Chessel, A. Dufour and J. Thioulouse. "The ade4 package-I-
One-table methods". In: _R News_ 4 (2004), pp. 5-10.

[4] A. Couture-Beil. _rjson: JSON for R_. R package version
0.2.15. 2014. <URL: https://CRAN.R-project.org/package=rjson>.

[5] S. Dray and A. Dufour. "The ade4 package: implementing the
duality diagram for ecologists". In: _Journal of Statistical
Software_ 22.4 (2007), pp. 1-20.

[6] S. Dray, A. Dufour and D. Chessel. "The ade4 package-II:
Two-table and K-table methods." In: _R News_ 7.2 (2007), pp.
47-52.

[7] L. Lahti, L. Elo, T. Aittokallio, et al. "Analysis of Probe
Reliability in Differential Gene Expression Studies with Short
Oligonucleotide Arrays". In: _TCBB/IEEE_ 8 (1 2011). R/BioC:
http://bioconductor.org/packages/release/bioc/html/RPA.html, pp.
217-225. <URL:
http://www.computer.org/portal/web/csdl/doi/10.1109/TCBB.2009.38.>.

[8] L. Lahti, A. Torrente, L. Elo, et al. "A fully scalable
online-preprocessing algorithm for short oligonucleotide
microarray atlases". In: _Nucleic Acids Research_ 41 (10 2013), p.
e110. <URL: URL:
http://nar.oxfordjournals.org/content/41/10/e110>.

[9] P. Langfelder and S. Horvath. "Fast R Functions for Robust
Correlations and Hierarchical Clustering". In: _Journal of
Statistical Software_ 46.11 (2012), pp. 1-17. <URL:
http://www.jstatsoft.org/v46/i11/>.

[10] P. Langfelder and S. Horvath. "WGCNA: an R package for
weighted correlation network analysis". In: _BMC Bioinformatics_
(2008), p. 559.

[11] P. J. McMurdie and S. Holmes. "phyloseq: An R package for
reproducible interactive analysis and graphics of microbiome
census data". In: _PLoS ONE_ 8.4 (2013), p. e61217. <URL:
http://dx.plos.org/10.1371/journal.pone.0061217>.

[12] P. E. Meyer, F. Lafitte and G. Bontempi. "MINET: An open
source R/Bioconductor Package for Mutual Information based Network
Inference". In: _BMC Bioinformatics_ 9 (2008). <URL:
http://www.biomedcentral.com/1471-2105/9/461>.

[13] D. Müllner. "fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python". In: _Journal of Statistical
Software_ 53.9 (2013), pp. 1-18. <URL:
http://www.jstatsoft.org/v53/i09/>.

[14] J. Oksanen, F. G. Blanchet, R. Kindt, et al. _vegan:
Community Ecology Package_. R package version 2.3-5. 2016. <URL:
https://CRAN.R-project.org/package=vegan>.

[15] R Core Team. _R: A language and environment for statistical
computing_. Vienna, Austria: R Foundation for Statistical
Computing, 2013. ISBN: ISBN 3-900051-07-0. <URL:
http://www.R-project.org/>.

[16] D. Temple Lang and t. C. team. _RCurl: General Network
(HTTP/FTP/...) Client Interface for R_. R package version
1.95-4.8. 2016. <URL: https://CRAN.R-project.org/package=RCurl>.

[17] H. Wickham. _ggplot2: Elegant Graphics for Data Analysis_.
Springer-Verlag New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://ggplot2.org>.

[18] Wickham and Hadley. "Reshaping data with the reshape
package". In: _Journal of Statistical Software_ 21.12 (2007).
<URL: http://www.jstatsoft.org/v21/i12/paper>.

[19] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.4.3. 2015. <URL:
https://CRAN.R-project.org/package=dplyr>.

### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=en_US.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=en_US.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=en_US.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=en_US.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] tcltk     parallel  splines   stats4    grid      stats     graphics 
##  [8] grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] HITChipDB_0.5.29      RPA_1.27.41           affy_1.48.0          
##  [4] Biobase_2.30.0        BiocGenerics_0.16.1   RMySQL_0.10.8        
##  [7] preprocessCore_1.32.0 DBI_0.3.1             FD_1.0-12            
## [10] geometry_0.3-6        magic_1.5-6           abind_1.4-3          
## [13] ape_3.4               ade4_1.7-4            microbiome_0.99.81   
## [16] vegan_2.3-5           permute_0.9-0         earlywarnings_1.1.22 
## [19] tseries_0.10-34       tgp_2.4-14            moments_0.14         
## [22] gridExtra_2.2.1       scales_0.4.0          SpiecEasi_0.1        
## [25] VGAM_1.0-0            huge_1.2.7            igraph_1.0.1         
## [28] lattice_0.20-33       Matrix_1.2-4          knitcitations_1.0.7  
## [31] knitr_1.12.3          intergraph_2.0-2      sna_2.3-2            
## [34] network_1.13.0        ggnet_0.1.0           GGally_1.0.1         
## [37] devtools_1.10.0       limma_3.26.9          sorvi_0.7.43         
## [40] tibble_1.0            ggplot2_2.1.0         tidyr_0.4.1          
## [43] dplyr_0.4.3           MASS_7.3-45           netresponse_1.20.15  
## [46] reshape2_1.4.1        mclust_5.2            minet_3.28.0         
## [49] Rgraphviz_2.14.0      graph_1.48.0          phyloseq_1.14.0      
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-126         bitops_1.0-6         lubridate_1.5.0     
##  [4] RColorBrewer_1.1-2   httr_1.1.0           tools_3.2.3         
##  [7] affyio_1.40.0        R6_2.1.2             rpart_4.1-10        
## [10] KernSmooth_2.23-15   dmt_0.8.20           lazyeval_0.1.10     
## [13] nortest_1.0-4        mgcv_1.8-12          colorspace_1.2-6    
## [16] withr_1.0.1          chron_2.3-47         formatR_1.3         
## [19] labeling_0.3         lmtest_0.9-34        mvtnorm_1.0-5       
## [22] quadprog_1.5-5       stringr_1.0.0        digest_0.6.9        
## [25] XVector_0.10.0       bibtex_0.4.0         highr_0.5.1         
## [28] maps_3.1.0           BiocInstaller_1.20.1 zoo_1.7-12          
## [31] RCurl_1.95-4.8       magrittr_1.5         Rcpp_0.12.4         
## [34] munsell_0.4.3        S4Vectors_0.8.11     maptree_1.4-7       
## [37] RefManageR_0.10.6    stringi_1.0-1        RJSONIO_1.3-0       
## [40] zlibbioc_1.16.0      plyr_1.8.3           qvalue_2.2.2        
## [43] Biostrings_2.38.4    multtest_2.26.0      boot_1.3-18         
## [46] codetools_0.2-14     XML_3.98-1.4         evaluate_0.8.3      
## [49] biom_0.3.12          data.table_1.9.6     spam_1.3-0          
## [52] foreach_1.4.3        gtable_0.2.0         reshape_0.8.5       
## [55] assertthat_0.1       roxygen2_5.0.1       Kendall_2.2         
## [58] survival_2.38-3      iterators_1.0.8      som_0.3-5           
## [61] memoise_1.0.0        IRanges_2.4.8        fields_8.3-6        
## [64] cluster_2.0.3
```




