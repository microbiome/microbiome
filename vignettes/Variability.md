### Inter-individual variability (variability within group of samples)

Assess 'inter-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between samples and their mean for a given samples vs phylotypes matrix. For the illustration, calculate inter-individual stability (variability) separately for Placebo and LGG groups.

Load example data


```r
library(microbiome)
x <- download_microbiome("dietswap")
# Add time field (two time points needed within each group for the 
# intraindividual method)
sample_data(x)$time <- sample_data(x)$timepoint.within.group
```


Variability across subjects within a group


```r
res <- estimate_variability(x, "interindividual")
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Inter-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

![plot of chunk variability-example2d](figure/variability-example2d-1.png) 


### Intra-individual variability

Variability within subjects over time (intra-individual stability). Assess 'intra-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between two time points within subjects, calculated separately within each group. For illustration, check intra-individual stability (variability) separately for Placebo and LGG groups.


```r
res <- estimate_variability(x, "intraindividual")
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Intra-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

```
## Warning: Removed 6 rows containing non-finite values (stat_boxplot).
```

![plot of chunk variability-intra](figure/variability-intra-1.png) 


### Version information


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

