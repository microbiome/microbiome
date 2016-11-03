---
title: "Variability"
author: "Leo Lahti"
date: "2016-11-03"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - variability}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


### Inter-individual homogeneity (within group of samples)

Assess 'inter-individual stability', or homogeneity, as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between the samples and their mean for a given samples vs phylotypes matrix. For illustration, calculate inter-individual homogeneity separately for Placebo and LGG groups. Note that this homogeneity measure is affected by sample size.

Load example data


```r
library(microbiome)
data("dietswap")
x <- dietswap

# Add time field (two time points needed within each group for the 
# intraindividual method)
sample_data(x)$time <- sample_data(x)$timepoint.within.group
```


Heterogeneity across subjects within a group


```r
res <- estimate_homogeneity(x, "interindividual")
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Inter-individual homogeneity (p=", round(res$p.value, 6), ")", sep = ""))
p <- p + ylab("Correlation")
print(p)
```

![plot of chunk homogeneity-example2d](figure/homogeneity-example2d-1.png)


### Intra-individual stability

Homogeneity within subjects over time (also called intra-individual stability in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html)). Defined as the average correlation between two time points within subjects within each group. For illustration, check intra-individual stability (homogeneity) separately for Placebo and LGG groups.


```r
res <- estimate_homogeneity(x, "intraindividual")
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Intra-individual homogeneity (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

![plot of chunk homogeneity-intra](figure/homogeneity-intra-1.png)


### Time series


```r
data("atlas1006")
pseq <- atlas1006
pseq <- subset_samples(pseq, DNA_extraction_method == "r")
pseq <- transform_phyloseq(pseq, "relative.abundance")
p <- plot_timeseries(pseq, "Dialister", subject = "831", tipping.point = 0.5)
print(p)
```

![plot of chunk homogeneity-timeseries](figure/homogeneity-timeseries-1.png)

Pick samples at the baseline time points only:


```r
data("atlas1006")
pseq0 <- pick_baseline(atlas1006)
```


## Further visualization tools

Draw regression curve with smoothed error bars based on
the [Visually-Weighted Regression](http://www.fight-entropy.com/2012/07/visually-weighted-regression.html) by Solomon M. Hsiang. The sorvi implementation extends [Felix Schonbrodt's original code](http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/).


```r
data(atlas1006)
p <- plot_regression(diversity ~ age, sample_data(atlas1006))
print(p)
```

![plot of chunk variability-regression](figure/variability-regression-1.png)

### Version information


```r
sessionInfo()
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.1 LTS
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
##  [1] FD_1.0-12             vegan_2.4-1           lattice_0.20-34      
##  [4] permute_0.9-4         geometry_0.3-6        magic_1.5-6          
##  [7] abind_1.4-5           ape_3.5               ade4_1.7-4           
## [10] RColorBrewer_1.1-2    knitcitations_1.0.7.1 knitr_1.14           
## [13] intergraph_2.0-2      sna_2.4               statnet.common_3.3.0 
## [16] network_1.13.0        ggnet_0.1.0           GGally_1.2.0         
## [19] devtools_1.12.0       limma_3.28.21         sorvi_0.7.46         
## [22] tibble_1.2            ggplot2_2.1.0         tidyr_0.6.0          
## [25] dplyr_0.5.0           MASS_7.3-45           netresponse_1.32.2   
## [28] reshape2_1.4.1        mclust_5.2            minet_3.30.0         
## [31] Rgraphviz_2.16.0      graph_1.50.0          BiocGenerics_0.18.0  
## [34] microbiome_0.99.87    phyloseq_1.16.2      
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-128          bitops_1.0-6          matrixStats_0.51.0   
##  [4] lubridate_1.6.0       httr_1.2.1            doParallel_1.0.10    
##  [7] dynamicTreeCut_1.63-1 tools_3.3.2           R6_2.2.0             
## [10] rpart_4.1-10          dmt_0.8.20            lazyeval_0.2.0       
## [13] Hmisc_3.17-4          DBI_0.5-1             mgcv_1.8-15          
## [16] colorspace_1.2-7      nnet_7.3-12           withr_1.0.2          
## [19] gridExtra_2.2.1       moments_0.14          preprocessCore_1.34.0
## [22] chron_2.3-47          WGCNA_1.51            Biobase_2.32.0       
## [25] formatR_1.4           Cairo_1.5-9           labeling_0.3         
## [28] scales_0.4.0          mvtnorm_1.0-5         tgp_2.4-14           
## [31] stringr_1.1.0         digest_0.6.10         foreign_0.8-67       
## [34] XVector_0.12.1        bibtex_0.4.0          RSQLite_1.0.0        
## [37] impute_1.46.0         jsonlite_1.1          acepack_1.3-3.3      
## [40] RCurl_1.95-4.8        magrittr_1.5          GO.db_3.3.0          
## [43] Formula_1.2-1         biomformat_1.0.2      Matrix_1.2-7.1       
## [46] Rcpp_0.12.7           munsell_0.4.3         S4Vectors_0.10.3     
## [49] maptree_1.4-7         RefManageR_0.11.0     stringi_1.1.2        
## [52] RJSONIO_1.3-0         zlibbioc_1.18.0       rhdf5_2.16.0         
## [55] plyr_1.8.4            qvalue_2.4.2          Biostrings_2.40.2    
## [58] splines_3.3.2         multtest_2.28.0       igraph_1.0.1         
## [61] fastcluster_1.1.21    codetools_0.2-15      stats4_3.3.2         
## [64] XML_3.98-1.4          evaluate_0.10         latticeExtra_0.6-28  
## [67] data.table_1.9.6      foreach_1.4.3         gtable_0.2.0         
## [70] reshape_0.8.5         assertthat_0.1        survival_2.39-5      
## [73] iterators_1.0.8       AnnotationDbi_1.34.4  memoise_1.0.0        
## [76] IRanges_2.6.1         cluster_2.0.5
```

