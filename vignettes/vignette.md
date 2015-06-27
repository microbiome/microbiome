---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2015-06-27"
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

For probe-level summarization of phylogenetic microarray data we use and recommend Robust Probabilistic Averaging (RPA) [(Lahti, Torrente, Elo, Brazma, and Rung, 2013; Lahti, Elo, Aittokallio, and Kaski, 2011); (Lahti, Elo, Aittokallio, and Kaski, 2011)].

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
 * mixOmics (Cao, Gonzalez, Rohart, Gautier, Monget, Coquery, Yao, and Liquet., 2015)
 * RCurl (Temple Lang, 2015)
 * vegan (Oksanen, Blanchet, Kindt, Legendre, Minchin, O'Hara, Simpson, Solymos, Stevens, and Wagner, 2015)
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
5.0-4. 2015. <URL: http://CRAN.R-project.org/package=mixOmics>.

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
Community Ecology Package_. R package version 2.3-0. 2015. <URL:
http://CRAN.R-project.org/package=vegan>.

[15] R Core Team. _R: A language and environment for statistical
computing_. Vienna, Austria: R Foundation for Statistical
Computing, 2013. ISBN: ISBN 3-900051-07-0. <URL:
http://www.R-project.org/>.

[16] D. Temple Lang. _RCurl: General network (HTTP/FTP/...) client
interface for R_. R package version 1.95-4.6. 2015. <URL:
http://CRAN.R-project.org/package=RCurl>.

[17] H. Wickham. _ggplot2: elegant graphics for data analysis_.
Springer New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://had.co.nz/ggplot2/book>.

[18] H. Wickham and R. Francois. _dplyr: A Grammar of Data
Manipulation_. R package version 0.4.2. 2015. <URL:
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
## R version 3.2.1 (2015-06-18)
## Platform: x86_64-unknown-linux-gnu (64-bit)
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
##  [1] gridExtra_0.9.1     googleVis_0.5.8     limma_3.24.10      
##  [4] vegan_2.3-0         lattice_0.20-31     permute_0.8-4      
##  [7] scales_0.2.5        RSQLite_1.0.0       DBI_0.3.1          
## [10] mgcv_1.8-6          nlme_3.1-120        dplyr_0.4.2        
## [13] netresponse_1.18.0  reshape_0.8.5       mclust_5.0.1       
## [16] minet_3.26.0        Rgraphviz_2.12.0    graph_1.46.0       
## [19] ggplot2_1.0.1       sorvi_0.7.26        microbiome_0.99.55 
## [22] RPA_1.24.0          affy_1.46.1         Biobase_2.28.0     
## [25] BiocGenerics_0.14.0 phyloseq_1.13.2     rdryad_0.1.1       
## [28] knitcitations_1.0.6 knitr_1.10.5        rmarkdown_0.7      
## [31] scimapClient_0.2.1 
## 
## loaded via a namespace (and not attached):
##   [1] spam_1.0-1                Hmisc_3.16-0             
##   [3] plyr_1.8.3                igraph_0.7.1             
##   [5] lazyeval_0.1.10           splines_3.2.1            
##   [7] BiocParallel_1.2.5        GenomeInfoDb_1.4.1       
##   [9] digest_0.6.8              foreach_1.4.2            
##  [11] BiocInstaller_1.18.3      htmltools_0.2.6          
##  [13] GO.db_3.1.2               gdata_2.16.1             
##  [15] magrittr_1.5              memoise_0.2.1            
##  [17] cluster_2.0.2             doParallel_1.0.8         
##  [19] fastcluster_1.1.16        Biostrings_2.36.1        
##  [21] annotate_1.46.0           matrixStats_0.14.1       
##  [23] Kendall_2.2               tseries_0.10-34          
##  [25] colorspace_1.2-6          RCurl_1.95-4.6           
##  [27] RcppArmadillo_0.5.200.1.0 genefilter_1.50.0        
##  [29] impute_1.42.0             survival_2.38-2          
##  [31] zoo_1.7-12                iterators_1.0.7          
##  [33] ape_3.3                   gtable_0.1.2             
##  [35] zlibbioc_1.14.0           XVector_0.8.0            
##  [37] RGCCA_2.0                 tgp_2.4-11               
##  [39] maps_2.3-9                futile.options_1.0.0     
##  [41] pheatmap_1.0.2            mvtnorm_1.0-2            
##  [43] som_0.3-5                 bibtex_0.4.0             
##  [45] Rcpp_0.11.6               xtable_1.7-4             
##  [47] foreign_0.8-63            preprocessCore_1.30.0    
##  [49] Formula_1.2-1             stats4_3.2.1             
##  [51] httr_0.6.1                RColorBrewer_1.1-2       
##  [53] acepack_1.3-3.3           XML_3.98-1.2             
##  [55] nnet_7.3-9                locfit_1.5-9.1           
##  [57] RJSONIO_1.3-0             dynamicTreeCut_1.62      
##  [59] labeling_0.3              reshape2_1.4.1           
##  [61] AnnotationDbi_1.30.1      munsell_0.4.2            
##  [63] tools_3.2.1               moments_0.14             
##  [65] ade4_1.7-2                evaluate_0.7             
##  [67] stringr_1.0.0             dmt_0.8.20               
##  [69] maptree_1.4-7             RefManageR_0.8.63        
##  [71] rgl_0.95.1247             formatR_1.2              
##  [73] affyio_1.36.0             geneplotter_1.46.0       
##  [75] stringi_0.5-2             highr_0.5                
##  [77] futile.logger_1.4.1       fields_8.2-1             
##  [79] earlywarnings_1.1.21      Matrix_1.2-1             
##  [81] multtest_2.24.0           biom_0.3.12              
##  [83] lmtest_0.9-34             data.table_1.9.4         
##  [85] bitops_1.0-6              GenomicRanges_1.20.5     
##  [87] qvalue_2.0.0              R6_2.0.1                 
##  [89] latticeExtra_0.6-26       KernSmooth_2.23-14       
##  [91] IRanges_2.2.4             codetools_0.2-11         
##  [93] lambda.r_1.1.7            boot_1.3-16              
##  [95] MASS_7.3-41               gtools_3.5.0             
##  [97] assertthat_0.1            chron_2.3-46             
##  [99] proto_0.3-10              DESeq2_1.8.1             
## [101] OAIHarvester_0.1-7        nortest_1.0-3            
## [103] S4Vectors_0.6.0           diptest_0.75-7           
## [105] mixOmics_5.0-4            quadprog_1.5-5           
## [107] rpart_4.1-9               WGCNA_1.47               
## [109] lubridate_1.3.3
```




