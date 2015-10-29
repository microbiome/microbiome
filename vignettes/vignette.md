---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2015-10-29"
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

The microbiome package contains general-purpose tools for microarray-based analysis of microbiome profiling data sets in R (R Core Team, 2013). 

The analyses are based on the phyloseq class from the independent [phyloseq](https://github.com/joey711/phyloseq) R package. Convert your data into phyloseq format before testing these examples.

For probe-level summarization of phylogenetic microarray data we use and recommend Robust Probabilistic Averaging (RPA) [(Lahti, Torrente, Elo, et al., 2013; Lahti, Elo, Aittokallio, et al., 2011); (Lahti, Elo, Aittokallio, et al., 2011b)].

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

### Visualization and related tools

Download some [example data sets](Data.md) to try these tools:

* [Barplots](Composition.md)
* [Boxplots](Boxplots.md)
* [Density](Density.md)
* [Heatmaps](Heatmap.md)
* [Motion charts](Motionchart.md)
* [Networks](Networks.md)
* [Ordination](Ordination.md)
* [PCA](Ordination.md)
* [Cross hybridization](Crosshyb.md)


### Clustering 
* [Bimodality](Stability.md)
* [Clustering](Clustering.md)


### Microbiota composition
* [Bistability analysis](Stability.md)
* [Composition](Composition.md)
* [Core microbiota](Core.md)
* [Diversity](Diversity.md)
* [Probe level studies](Probelevel.md)
* [Stability](Stability.md)
* [Tipping elements](Stability.md)
* [Variability](Variability.md) (Intra- and inter-individual 'stability')


### Linear models, comparisons, and association studies
* [Linear models](limma.md)
* [Pairwise comparisons](Comparisons.md)


### Other statistical analysis
* [ROC curves](ROC.md)
* [RDA](RDA.md)


### Output
* [Generating output files/figures](Output.md)


### Misc
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
 * mixOmics (Cao, Gonzalez, Rohart, et al., 2015)
 * RCurl (Temple Lang and team, 2015)
 * vegan (Oksanen, Blanchet, Kindt, et al., 2015)
 * reshape (Wickham and Hadley, 2007)
 * WGCNA (Langfelder and Horvath, 2008; Langfelder and Horvath, 2012)
 * ggplot2 (Wickham, 2009)
 * RPA (Lahti, Torrente, Elo, et al., 2013; Lahti, Elo, Aittokallio, et al., 2011)
 * minet (Meyer, Lafitte, and Bontempi, 2008)
 * fastcluster (Müllner, 2013)
 * dplyr (Wickham and Francois, 2015)


### References



[1] N. Caballero. _df2json: Convert a dataframe to JSON_. R
package version 0.0.2. 2013. <URL:
http://CRAN.R-project.org/package=df2json>.

[2] K. L. Cao, I. Gonzalez, S. D. w. k. c. F. Rohart, et al.
_mixOmics: Omics Data Integration Project_. R package version
5.1.2. 2015. <URL: http://CRAN.R-project.org/package=mixOmics>.

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

[7] L. Lahti, L. L. Elo, T. Aittokallio, et al. "Probabilistic
Analysis of Probe Reliability in Differential Gene Expression
Studies with Short Oligonucleotide Arrays". In: _IEEE/ACM Trans.
Comput. Biol. and Bioinf._ 8.1 (Jan. 2011), pp. 217-225. DOI:
10.1109/tcbb.2009.38. <URL:
http://dx.doi.org/10.1109/TCBB.2009.38>.

[8] L. Lahti, L. Elo, T. Aittokallio, et al. "Analysis of Probe
Reliability in Differential Gene Expression Studies with Short
Oligonucleotide Arrays". In: _TCBB/IEEE_ 8 (1 2011). R/BioC:
http://bioconductor.org/packages/release/bioc/html/RPA.html, pp.
217-225. <URL:
http://www.computer.org/portal/web/csdl/doi/10.1109/TCBB.2009.38.>.

[9] L. Lahti, A. Torrente, L. Elo, et al. "A fully scalable
online-preprocessing algorithm for short oligonucleotide
microarray atlases". In: _Nucleic Acids Research_ 41 (10 2013), p.
e110. <URL: URL:
http://nar.oxfordjournals.org/content/41/10/e110>.

[10] P. Langfelder and S. Horvath. "Fast R Functions for Robust
Correlations and Hierarchical Clustering". In: _Journal of
Statistical Software_ 46.11 (2012), pp. 1-17. <URL:
http://www.jstatsoft.org/v46/i11/>.

[11] P. Langfelder and S. Horvath. "WGCNA: an R package for
weighted correlation network analysis". In: _BMC Bioinformatics_
(2008), p. 559.

[12] P. E. Meyer, F. Lafitte and G. Bontempi. "MINET: An open
source R/Bioconductor Package for Mutual Information based Network
Inference". In: _BMC Bioinformatics_ 9 (2008). <URL:
http://www.biomedcentral.com/1471-2105/9/461>.

[13] D. Müllner. "fastcluster: Fast Hierarchical, Agglomerative
Clustering Routines for R and Python". In: _Journal of Statistical
Software_ 53.9 (2013), pp. 1-18. <URL:
http://www.jstatsoft.org/v53/i09/>.

[14] J. Oksanen, F. G. Blanchet, R. Kindt, et al. _vegan:
Community Ecology Package_. R package version 2.3-1. 2015. <URL:
http://CRAN.R-project.org/package=vegan>.

[15] R Core Team. _R: A language and environment for statistical
computing_. Vienna, Austria: R Foundation for Statistical
Computing, 2013. ISBN: ISBN 3-900051-07-0. <URL:
http://www.R-project.org/>.

[16] D. Temple Lang and t. C. team. _RCurl: General Network
(HTTP/FTP/...) Client Interface for R_. R package version
1.95-4.7. 2015. <URL: http://CRAN.R-project.org/package=RCurl>.

[17] H. Wickham. _ggplot2: elegant graphics for data analysis_.
Springer New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://had.co.nz/ggplot2/book>.

[18] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.4.3. 2015. <URL:
http://CRAN.R-project.org/package=dplyr>.

[19] Wickham and Hadley. "Reshaping data with the reshape
package". In: _Journal of Statistical Software_ 21.12 (2007).
<URL: http://www.jstatsoft.org/v21/i12/paper>.

### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 15.04
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
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] googleVis_0.5.10    gridExtra_2.0.0     limma_3.24.15      
##  [4] mgcv_1.8-8          nlme_3.1-122        dplyr_0.4.3        
##  [7] netresponse_1.18.0  reshape_0.8.5       mclust_5.1         
## [10] minet_3.26.0        Rgraphviz_2.12.0    graph_1.46.0       
## [13] ggplot2_1.0.1       sorvi_0.7.32        microbiome_0.99.63 
## [16] RPA_1.24.0          affy_1.46.1         Biobase_2.28.0     
## [19] BiocGenerics_0.14.0 phyloseq_1.12.2     rdryad_0.1.1       
## [22] knitcitations_1.0.7 knitr_1.11         
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.2-6          qvalue_2.0.0             
##   [3] som_0.3-5                 futile.logger_1.4.1      
##   [5] XVector_0.8.0             OAIHarvester_0.1-7       
##   [7] RcppArmadillo_0.6.100.0.0 GenomicRanges_1.20.8     
##   [9] affyio_1.36.0             mvtnorm_1.0-3            
##  [11] AnnotationDbi_1.30.1      lubridate_1.3.3          
##  [13] RefManageR_0.8.63         codetools_0.2-14         
##  [15] splines_3.2.2             geneplotter_1.46.0       
##  [17] mixOmics_5.1.2            tgp_2.4-11               
##  [19] ade4_1.7-2                spam_1.3-0               
##  [21] Formula_1.2-1             annotate_1.46.1          
##  [23] cluster_2.0.3             pheatmap_1.0.7           
##  [25] Kendall_2.2               httr_1.0.0               
##  [27] lazyeval_0.1.10           assertthat_0.1           
##  [29] Matrix_1.2-2              formatR_1.2.1            
##  [31] acepack_1.3-3.3           tools_3.2.2              
##  [33] igraph_1.0.1              gtable_0.1.2             
##  [35] reshape2_1.4.1            maps_3.0.0-2             
##  [37] Rcpp_0.12.1               Biostrings_2.36.4        
##  [39] RJSONIO_1.3-0             multtest_2.24.0          
##  [41] biom_0.3.12               gdata_2.17.0             
##  [43] ape_3.3                   preprocessCore_1.30.0    
##  [45] iterators_1.0.8           lmtest_0.9-34            
##  [47] stringr_1.0.0             proto_0.3-10             
##  [49] gtools_3.5.0              XML_3.98-1.3             
##  [51] zlibbioc_1.14.0           MASS_7.3-44              
##  [53] zoo_1.7-12                scales_0.3.0             
##  [55] BiocInstaller_1.18.5      lambda.r_1.1.7           
##  [57] RColorBrewer_1.1-2        fields_8.3-5             
##  [59] memoise_0.2.1             rpart_4.1-10             
##  [61] latticeExtra_0.6-26       stringi_1.0-1            
##  [63] maptree_1.4-7             RSQLite_1.0.0            
##  [65] highr_0.5.1               genefilter_1.50.0        
##  [67] S4Vectors_0.6.6           tseries_0.10-34          
##  [69] foreach_1.4.3             nortest_1.0-4            
##  [71] permute_0.8-4             boot_1.3-17              
##  [73] BiocParallel_1.2.22       bibtex_0.4.0             
##  [75] chron_2.3-47              GenomeInfoDb_1.4.3       
##  [77] moments_0.14              bitops_1.0-6             
##  [79] dmt_0.8.20                rgl_0.95.1367            
##  [81] evaluate_0.8              lattice_0.20-33          
##  [83] labeling_0.3              plyr_1.8.3               
##  [85] magrittr_1.5              DESeq2_1.8.2             
##  [87] R6_2.1.1                  IRanges_2.2.9            
##  [89] earlywarnings_1.1.22      Hmisc_3.17-0             
##  [91] DBI_0.3.1                 foreign_0.8-66           
##  [93] survival_2.38-3           RCurl_1.95-4.7           
##  [95] nnet_7.3-11               futile.options_1.0.0     
##  [97] KernSmooth_2.23-15        ellipse_0.3-8            
##  [99] locfit_1.5-9.1            data.table_1.9.6         
## [101] vegan_2.3-1               digest_0.6.8             
## [103] diptest_0.75-7            xtable_1.7-4             
## [105] stats4_3.2.2              munsell_0.4.2            
## [107] quadprog_1.5-5
```




