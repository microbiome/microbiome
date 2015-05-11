---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2015-05-11"
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

The microbiome package contains general-purpose tools for microarray-based analysis of microbiome profiling data sets in R (R Core Team, 2013). See also the related but independent R package [phyloseq](https://github.com/joey711/phyloseq).


### How to start 

Example on reproducible document generation:  

 * [Minimal example](Template.md)  
 * [Example workflow](Atlas.md)

### Installation, example data sets and preprocessing
* [Installing R tools](Installation.md)

### Data
* [Download example data sets](Data.md)
* [Extract data from HITChip database](https://github.com/microbiome/HITChipDB/blob/master/vignettes/vignette.md)
* [Get phylogeny](Phylogeny.md)

### Visualization and related tools

Download some [example data sets](Data.md) to try these tools:

* [Barplots](Barplots.md)
* [Boxplots](Boxplots.md)
* [Density](Density.md)
* [Heatmaps](Heatmap.md)
* [Motion charts](Motionchart.md)
* [Oligo heatmap](Oligoheatmap.md)
* [PCA and other ordination/projection methods](PCA.md)
* [Cross hybridization](Crosshyb.md)

### Clustering 
* [Bimodality](Bimodality.md)
* [Clustering](Clustering.md)

### Microbiota composition
* [Core microbiota](Core.md)
* [Diversity](Diversity.md)
* [Probe level studies](Probelevel.md)
* [RelativeAbundance](RelativeAbundance.md)
* [Variability](Variability.md)
* [Bistability analysis](Stability.md)

### Linear models, comparisons, and association studies
* [Linear models](limma.md)
* [Pairwise comparisons](Comparisons.md)
* [Cross correlations](Crosscorrelation.md)

### Other statistical analysis
* [ROC curves](ROC.md)
* [RDA](RDA.md)
* [Phyloseq tools](Phyloseq.md)

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
 * mixOmics (Dejean, Gonzalez, Monget, Coquery, Yao, Liquet, and Rohart, 2014)
 * RCurl (Temple Lang, 2014)
 * vegan (Oksanen, Blanchet, Kindt, Legendre, Minchin, O'Hara, Simpson, Solymos, Stevens, and Wagner, 2015)
 * reshape (Wickham and Hadley, 2007)
 * WGCNA (Langfelder and Horvath, 2008; Langfelder and Horvath, 2012)
 * ggplot2 (Wickham, 2009)
 * RPA (Lahti, Torrente, Elo, et al., 2013; Lahti, Elo, Aittokallio, et al., 2011)
 * minet (Meyer, Lafitte, and Bontempi, 2008)
 * fastcluster (Müllner, 2013)
 * plyr (Wickham, 2011)


### References



[1] N. Caballero. _df2json: Convert a dataframe to JSON_. R
package version 0.0.2. 2013. <URL:
http://CRAN.R-project.org/package=df2json>.

[2] D. Chessel, A. Dufour and J. Thioulouse. "The ade4 package-I-
One-table methods". In: _R News_ 4 (2004), pp. 5-10.

[3] A. Couture-Beil. _rjson: JSON for R_. R package version
0.2.15. 2014. <URL: http://CRAN.R-project.org/package=rjson>.

[4] S. Dejean, I. Gonzalez, K. L. C. w. c. f. P. Monget, et al.
_mixOmics: Omics Data Integration Project_. R package version
5.0-3. 2014. <URL: http://CRAN.R-project.org/package=mixOmics>.

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
Community Ecology Package_. R package version 2.2-1. 2015. <URL:
http://CRAN.R-project.org/package=vegan>.

[14] R Core Team. _R: A language and environment for statistical
computing_. Vienna, Austria: R Foundation for Statistical
Computing, 2013. ISBN: ISBN 3-900051-07-0. <URL:
http://www.R-project.org/>.

[15] D. Temple Lang. _RCurl: General network (HTTP/FTP/...) client
interface for R_. R package version 1.95-4.5. 2014. <URL:
http://CRAN.R-project.org/package=RCurl>.

[16] H. Wickham. _ggplot2: elegant graphics for data analysis_.
Springer New York, 2009. ISBN: 978-0-387-98140-6. <URL:
http://had.co.nz/ggplot2/book>.

[17] H. Wickham. "The Split-Apply-Combine Strategy for Data
Analysis". In: _Journal of Statistical Software_ 40.1 (2011), pp.
1-29. <URL: http://www.jstatsoft.org/v40/i01/>.

[18] Wickham and Hadley. "Reshaping data with the reshape
package". In: _Journal of Statistical Software_ 21.12 (2007).
<URL: http://www.jstatsoft.org/v21/i12/paper>.

### Session info

This vignette was created with


```r
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
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
##  [1] tcltk     parallel  grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] scales_0.2.4          HITChipDB_0.5.15      RPA_1.24.0           
##  [4] affy_1.46.0           Biobase_2.28.0        BiocGenerics_0.14.0  
##  [7] RMySQL_0.10.3         DBI_0.3.1             preprocessCore_1.30.0
## [10] limma_3.24.0          gdata_2.13.3          netresponse_1.18.0   
## [13] reshape_0.8.5         mclust_5.0.0          minet_3.26.0         
## [16] Rgraphviz_2.12.0      graph_1.46.0          ggplot2_1.0.1        
## [19] microbiome_0.99.46    vegan_2.2-1           lattice_0.20-31      
## [22] permute_0.8-3         reshape2_1.4.1        phyloseq_1.11.2      
## [25] e1071_1.6-4           dplyr_0.4.1           ade4_1.7-2           
## [28] rdryad_0.1.1          knitcitations_1.0.5   knitr_1.9            
## [31] rmarkdown_0.5.1      
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.2-6        rjson_0.2.15           
##   [3] class_7.3-12            qvalue_2.0.0           
##   [5] som_0.3-5               futile.logger_1.4.1    
##   [7] XVector_0.8.0           OAIHarvester_0.1-7     
##   [9] RcppArmadillo_0.5.000.0 GenomicRanges_1.20.2   
##  [11] affyio_1.36.0           mvtnorm_1.0-2          
##  [13] AnnotationDbi_1.30.0    lubridate_1.3.3        
##  [15] RefManageR_0.8.45       codetools_0.2-11       
##  [17] splines_3.2.0           geneplotter_1.46.0     
##  [19] mixOmics_5.0-3          tgp_2.4-11             
##  [21] spam_1.0-1              Formula_1.2-1          
##  [23] annotate_1.46.0         cluster_2.0.1          
##  [25] pheatmap_1.0.2          Kendall_2.2            
##  [27] httr_0.6.1              lazyeval_0.1.10        
##  [29] assertthat_0.1          Matrix_1.2-0           
##  [31] formatR_1.2             acepack_1.3-3.3        
##  [33] htmltools_0.2.6         tools_3.2.0            
##  [35] igraph_0.7.1            gtable_0.1.2           
##  [37] maps_2.3-9              Rcpp_0.11.5            
##  [39] df2json_0.0.2           Biostrings_2.36.0      
##  [41] RJSONIO_1.3-0           multtest_2.24.0        
##  [43] biom_0.3.12             ape_3.2                
##  [45] nlme_3.1-120            iterators_1.0.7        
##  [47] lmtest_0.9-33           fastcluster_1.1.16     
##  [49] stringr_0.6.2           proto_0.3-10           
##  [51] gtools_3.4.2            devtools_1.7.0         
##  [53] XML_3.98-1.1            zlibbioc_1.14.0        
##  [55] MASS_7.3-40             zoo_1.7-12             
##  [57] BiocInstaller_1.18.1    lambda.r_1.1.7         
##  [59] RColorBrewer_1.1-2      fields_8.2-1           
##  [61] memoise_0.2.1           rpart_4.1-9            
##  [63] latticeExtra_0.6-26     maptree_1.4-7          
##  [65] RSQLite_1.0.0           genefilter_1.50.0      
##  [67] S4Vectors_0.6.0         tseries_0.10-34        
##  [69] foreach_1.4.2           nortest_1.0-3          
##  [71] boot_1.3-16             BiocParallel_1.2.0     
##  [73] bibtex_0.4.0            chron_2.3-45           
##  [75] GenomeInfoDb_1.4.0      moments_0.14           
##  [77] bitops_1.0-6            dmt_0.8.20             
##  [79] rgl_0.95.1201           evaluate_0.7           
##  [81] labeling_0.3            plyr_1.8.2             
##  [83] magrittr_1.5            DESeq2_1.8.0           
##  [85] IRanges_2.2.1           earlywarnings_1.1.19   
##  [87] Hmisc_3.15-0            foreign_0.8-63         
##  [89] mgcv_1.8-6              survival_2.38-1        
##  [91] RCurl_1.95-4.5          nnet_7.3-9             
##  [93] futile.options_1.0.0    KernSmooth_2.23-14     
##  [95] RGCCA_2.0               locfit_1.5-9.1         
##  [97] data.table_1.9.4        digest_0.6.8           
##  [99] xtable_1.7-4            stats4_3.2.0           
## [101] munsell_0.4.2           quadprog_1.5-5
```




