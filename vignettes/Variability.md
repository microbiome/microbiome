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
##  [4] RSQLite_1.0.0       DBI_0.3.1           mgcv_1.8-6         
##  [7] nlme_3.1-120        dplyr_0.4.2         netresponse_1.18.0 
## [10] reshape_0.8.5       mclust_5.0.1        minet_3.26.0       
## [13] Rgraphviz_2.12.0    graph_1.46.0        ggplot2_1.0.1      
## [16] sorvi_0.7.26        microbiome_0.99.60  RPA_1.24.0         
## [19] affy_1.46.1         Biobase_2.28.0      BiocGenerics_0.14.0
## [22] phyloseq_1.13.2     rdryad_0.1.1        knitcitations_1.0.6
## [25] knitr_1.10.5        rmarkdown_0.7       scimapClient_0.2.1 
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
##  [39] maps_2.3-9                scales_0.2.5             
##  [41] futile.options_1.0.0      pheatmap_1.0.2           
##  [43] mvtnorm_1.0-2             som_0.3-5                
##  [45] bibtex_0.4.0              Rcpp_0.11.6              
##  [47] xtable_1.7-4              foreign_0.8-63           
##  [49] preprocessCore_1.30.0     Formula_1.2-1            
##  [51] stats4_3.2.1              httr_0.6.1               
##  [53] RColorBrewer_1.1-2        acepack_1.3-3.3          
##  [55] XML_3.98-1.3              nnet_7.3-9               
##  [57] locfit_1.5-9.1            RJSONIO_1.3-0            
##  [59] dynamicTreeCut_1.62       labeling_0.3             
##  [61] reshape2_1.4.1            AnnotationDbi_1.30.1     
##  [63] munsell_0.4.2             tools_3.2.1              
##  [65] moments_0.14              ade4_1.7-2               
##  [67] evaluate_0.7              stringr_1.0.0            
##  [69] dmt_0.8.20                maptree_1.4-7            
##  [71] RefManageR_0.8.63         rgl_0.95.1247            
##  [73] formatR_1.2               affyio_1.36.0            
##  [75] geneplotter_1.46.0        stringi_0.5-5            
##  [77] highr_0.5                 futile.logger_1.4.1      
##  [79] fields_8.2-1              earlywarnings_1.1.22     
##  [81] lattice_0.20-31           Matrix_1.2-1             
##  [83] vegan_2.3-0               permute_0.8-4            
##  [85] multtest_2.24.0           biom_0.3.12              
##  [87] lmtest_0.9-34             data.table_1.9.4         
##  [89] bitops_1.0-6              GenomicRanges_1.20.5     
##  [91] qvalue_2.0.0              R6_2.0.1                 
##  [93] latticeExtra_0.6-26       KernSmooth_2.23-14       
##  [95] IRanges_2.2.4             codetools_0.2-11         
##  [97] lambda.r_1.1.7            boot_1.3-16              
##  [99] MASS_7.3-41               gtools_3.5.0             
## [101] assertthat_0.1            chron_2.3-47             
## [103] proto_0.3-10              DESeq2_1.8.1             
## [105] OAIHarvester_0.1-7        nortest_1.0-3            
## [107] S4Vectors_0.6.0           diptest_0.75-7           
## [109] mixOmics_5.0-4            quadprog_1.5-5           
## [111] rpart_4.1-9               WGCNA_1.47               
## [113] lubridate_1.3.3
```

