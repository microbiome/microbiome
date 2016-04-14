---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2016-04-13"
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


```
## Warning: All formats failed to parse. No formats found.
```

```
## Warning in ProcessDate(bib[["year"]], bib[["month"]]): Failed to parse
## month: July. Ignoring and using year only.
```


microbiome R package
===========

The microbiome package contains general-purpose tools for microarray-based analysis of microbiome profiling data sets in R (R Core Team, 2013). The are based on the phyloseq class from the independent [phyloseq](https://github.com/joey711/phyloseq) R package. The microbiome package adds extra functionality for microbiomics data sets that can be represented in this format.


### How to start 

Example on reproducible document generation:  

 * [Installing R tools](Installation.md)
 * [How to get started](Template.md)  
 * [Example workflow](Atlas.md)


### Data

* [Download example data sets](Data.md)
* [Extract data from HITChip database](https://github.com/microbiome/HITChipDB/blob/master/vignettes/vignette.md)
* [Taxonomy](Phylogeny.md)
* [Preprocessing and Filtering](Preprocessing.md)
* [Importing data from different formats](Import.md)

### Visualization and related tools

Download some [example data sets](Data.md) to try these tools:

* [Barplots](Composition.md)
* [Boxplots](Boxplots.md)
* [Density](Density.md)
* [Heatmaps](Heatmap.md)
* [Networks](Networks.md)

### Interactive  

* [Interactive](Interactive.md)
* [Motion charts](Motionchart.md)

### Ordination  

* [RDA](RDA.md)
* [Ordination](Ordination.md)
* [PCA](Ordination.md)

### Clustering 
* [Bimodality](Stability.md)
* [Clustering](Clustering.md)

### Microbiota composition
* [Bistability analysis](Stability.md)
* [Composition](Composition.md)
* [Core microbiota](Core.md)
* [Diversity](Diversity.md)
* [Indices](Indices.md)
* [Stability](Stability.md)
* [Tipping elements](Stability.md)
* [Variability](Variability.md) (Intra- and inter-individual 'stability')


### Linear models, comparisons, and association studies
* [Linear models](limma.md)
* [Pairwise comparisons](Comparisons.md)


### Phylogenetic microarrays  

* [Cross hybridization](Crosshyb.md)
* [Probe level studies](Probelevel.md)


### Misc

* [Generating output files/figures](Output.md)
* [Miscellaneous](misc.md)


### Licensing and Citations

This work can be freely used, modified and distributed under the 
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD\_licenses).

Kindly cite the work as 'Leo Lahti and Jarkko Salojarvi
(2014). microbiome R package. URL: http://microbiome.github.com'.


### Dependencies

The package utilizes tools from a number of other CRAN and
Bioconductor extensions, including:

 * df2json (Caballero, 2013)
 * rjson (Couture-Beil, 2014)
 * ade4 (Dray and Dufour, 2007; Chessel, Dufour, and Thioulouse, 2004; Dray, Dufour, and Chessel, 2007)
 * mixOmics (Cao, Gonzalez, Rohart, Gautier, Monget, Coquery, Yao, and Liquet., 2015)
 * RCurl (Temple Lang and team, 2016)
 * vegan (Oksanen, Blanchet, Kindt, Legendre, Minchin, O'Hara, Simpson, Solymos, Stevens, and Wagner, 2016)
 * reshape (Wickham and Hadley, 2007)
 * WGCNA (Langfelder and Horvath, 2008; Langfelder and Horvath, 2012)
 * ggplot2 (Wickham, 2009)
 * RPA (Lahti, Torrente, Elo, Brazma, and Rung, 2013; Lahti, Elo, Aittokallio, and Kaski, 2011)
 * minet (Meyer, Lafitte, and Bontempi, 2008)
 * fastcluster (Müllner, 2013)
 * dplyr (Wickham and Francois, 2015)


### References



[1] N. Caballero. _df2json: Convert a dataframe to JSON_. R
package version 0.0.2. 2013. <URL:
http://CRAN.R-project.org/package=df2json>.

[2] K. L. Cao, I. Gonzalez, S. D. w. k. c. F. Rohart, et al.
_mixOmics: Omics Data Integration Project_. R package version
5.2.0. 2015. <URL: http://CRAN.R-project.org/package=mixOmics>.

[3] D. Chessel, A. Dufour and J. Thioulouse. "The ade4 package-I-
One-table methods". In: _R News_ 4 (2004), pp. 5-10.

[4] A. Couture-Beil. _rjson: JSON for R_. R package version
0.2.15. 2014. <URL: http://CRAN.R-project.org/package=rjson>.

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

[11] P. E. Meyer, F. Lafitte and G. Bontempi. "MINET: An open
source R/Bioconductor Package for Mutual Information based Network
Inference". In: _BMC Bioinformatics_ 9 (2008). <URL:
http://www.biomedcentral.com/1471-2105/9/461>.

[12] D. Müllner. "fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python". In: _Journal of Statistical
Software_ 53.9 (2013), pp. 1-18. <URL:
http://www.jstatsoft.org/v53/i09/>.

[13] J. Oksanen, F. G. Blanchet, R. Kindt, et al. _vegan:
Community Ecology Package_. R package version 2.3-5. 2016. <URL:
http://CRAN.R-project.org/package=vegan>.

[14] R Core Team. _R: A language and environment for statistical
computing_. Vienna, Austria: R Foundation for Statistical
Computing, 2013. ISBN: ISBN 3-900051-07-0. <URL:
http://www.R-project.org/>.

[15] D. Temple Lang and t. C. team. _RCurl: General Network
(HTTP/FTP/...) Client Interface for R_. R package version
1.95-4.8. 2016. <URL: http://CRAN.R-project.org/package=RCurl>.

[16] H. Wickham. _ggplot2: Elegant Graphics for Data Analysis_.
Springer-Verlag New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://ggplot2.org>.

[17] Wickham and Hadley. "Reshaping data with the reshape
package". In: _Journal of Statistical Software_ 21.12 (2007).
<URL: http://www.jstatsoft.org/v21/i12/paper>.

[18] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.4.3. 2015. <URL:
http://CRAN.R-project.org/package=dplyr>.

### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 15.10
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
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] vegan_2.3-5          lattice_0.20-33      permute_0.9-0       
##  [4] RSQLite_1.0.0        DBI_0.3.1            earlywarnings_1.1.22
##  [7] tseries_0.10-34      tgp_2.4-14           moments_0.14        
## [10] gridExtra_2.2.1      googleVis_0.5.10     scales_0.4.0        
## [13] knitcitations_1.0.7  knitr_1.12.3         intergraph_2.0-2    
## [16] sna_2.3-2            network_1.13.0       ggnet_0.1.0         
## [19] GGally_1.0.1         devtools_1.11.0      limma_3.26.9        
## [22] sorvi_0.7.38         ggplot2_2.1.0        tidyr_0.4.1         
## [25] dplyr_0.4.3          MASS_7.3-45          netresponse_1.20.15 
## [28] reshape2_1.4.1       mclust_5.2           minet_3.28.0        
## [31] Rgraphviz_2.14.0     graph_1.48.0         microbiome_0.99.77  
## [34] RPA_1.27.3           phyloseq_1.14.0      affy_1.48.0         
## [37] Biobase_2.30.0       BiocGenerics_0.16.1 
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      rjson_0.2.15          dynamicTreeCut_1.63-1
##  [4] som_0.3-5             qvalue_2.2.2          XVector_0.10.0       
##  [7] affyio_1.40.0         AnnotationDbi_1.32.3  mvtnorm_1.0-5        
## [10] lubridate_1.5.6       RefManageR_0.10.13    xml2_0.1.2           
## [13] codetools_0.2-14      splines_3.2.2         doParallel_1.0.10    
## [16] impute_1.44.0         ade4_1.7-4            Formula_1.2-1        
## [19] spam_1.3-0            WGCNA_1.51            cluster_2.0.3        
## [22] GO.db_3.2.2           Kendall_2.2           oai_0.2.0            
## [25] httr_1.1.0            assertthat_0.1        Matrix_1.2-4         
## [28] lazyeval_0.1.10       formatR_1.3           acepack_1.3-3.3      
## [31] tools_3.2.2           igraph_1.0.1          rdryad_0.2.0         
## [34] gtable_0.2.0          maps_3.1.0            Rcpp_0.12.4          
## [37] Biostrings_2.38.4     RJSONIO_1.3-0         multtest_2.26.0      
## [40] biom_0.3.12           ape_3.4               preprocessCore_1.32.0
## [43] nlme_3.1-126          iterators_1.0.8       lmtest_0.9-34        
## [46] fastcluster_1.1.20    stringr_1.0.0         XML_3.98-1.4         
## [49] zlibbioc_1.16.0       zoo_1.7-12            BiocInstaller_1.20.1 
## [52] solr_0.1.6            RColorBrewer_1.1-2    fields_8.3-6         
## [55] memoise_1.0.0         rpart_4.1-10          latticeExtra_0.6-28  
## [58] reshape_0.8.5         stringi_1.0-1         maptree_1.4-7        
## [61] highr_0.5.1           S4Vectors_0.8.11      foreach_1.4.3        
## [64] nortest_1.0-4         boot_1.3-18           bibtex_0.4.0         
## [67] chron_2.3-47          bitops_1.0-6          matrixStats_0.50.1   
## [70] dmt_0.8.20            evaluate_0.8.3        labeling_0.3         
## [73] plyr_1.8.3            magrittr_1.5          R6_2.1.2             
## [76] IRanges_2.4.8         Hmisc_3.17-3          foreign_0.8-66       
## [79] withr_1.0.1           mgcv_1.8-12           nnet_7.3-12          
## [82] survival_2.38-3       RCurl_1.95-4.8        KernSmooth_2.23-15   
## [85] data.table_1.9.6      digest_0.6.9          diptest_0.75-7       
## [88] stats4_3.2.2          munsell_0.4.3         quadprog_1.5-5
```




