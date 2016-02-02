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


### Version information


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
##  [1] FD_1.0-12            geometry_0.3-6       magic_1.5-6         
##  [4] abind_1.4-3          ape_3.4              ade4_1.7-3          
##  [7] rdryad_0.2.0         RSQLite_1.0.0        DBI_0.3.1           
## [10] earlywarnings_1.1.22 tseries_0.10-34      tgp_2.4-11          
## [13] moments_0.14         gridExtra_2.0.0      vegan_2.3-3         
## [16] lattice_0.20-33      permute_0.8-4        knitcitations_1.0.7 
## [19] knitr_1.12           devtools_1.9.1       limma_3.26.5        
## [22] sorvi_0.7.35         ggplot2_2.0.0        tidyr_0.3.1         
## [25] dplyr_0.4.3          MASS_7.3-45          netresponse_1.21.14 
## [28] reshape2_1.4.1       mclust_5.1           minet_3.28.0        
## [31] Rgraphviz_2.14.0     graph_1.48.0         phyloseq_1.14.0     
## [34] microbiome_0.99.73   RPA_1.26.0           affy_1.48.0         
## [37] Biobase_2.30.0       BiocGenerics_0.16.1 
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      rjson_0.2.15          dynamicTreeCut_1.62  
##  [4] som_0.3-5             qvalue_2.2.2          XVector_0.10.0       
##  [7] affyio_1.40.0         AnnotationDbi_1.32.3  mvtnorm_1.0-3        
## [10] lubridate_1.5.0       RefManageR_0.10.5     xml2_0.1.2           
## [13] codetools_0.2-14      splines_3.2.2         doParallel_1.0.10    
## [16] impute_1.44.0         spam_1.3-0            Formula_1.2-1        
## [19] WGCNA_1.49            cluster_2.0.3         GO.db_3.2.2          
## [22] Kendall_2.2           oai_0.1.0             httr_1.0.0           
## [25] assertthat_0.1        Matrix_1.2-3          lazyeval_0.1.10      
## [28] formatR_1.2.1         acepack_1.3-3.3       tools_3.2.2          
## [31] igraph_1.0.1          gtable_0.1.2          maps_3.0.2           
## [34] Rcpp_0.12.3           Biostrings_2.38.3     RJSONIO_1.3-0        
## [37] multtest_2.26.0       biom_0.3.12           preprocessCore_1.32.0
## [40] nlme_3.1-122          iterators_1.0.8       lmtest_0.9-34        
## [43] fastcluster_1.1.16    stringr_1.0.0         XML_3.98-1.3         
## [46] zlibbioc_1.16.0       zoo_1.7-12            scales_0.3.0         
## [49] BiocInstaller_1.20.1  solr_0.1.6            RColorBrewer_1.1-2   
## [52] fields_8.3-6          curl_0.9.5            memoise_0.2.1        
## [55] rpart_4.1-10          latticeExtra_0.6-26   stringi_1.0-1        
## [58] maptree_1.4-7         highr_0.5.1           S4Vectors_0.8.7      
## [61] foreach_1.4.3         nortest_1.0-4         boot_1.3-17          
## [64] bibtex_0.4.0          chron_2.3-47          bitops_1.0-6         
## [67] matrixStats_0.50.1    dmt_0.8.20            evaluate_0.8         
## [70] labeling_0.3          plyr_1.8.3            magrittr_1.5         
## [73] R6_2.1.2              IRanges_2.4.6         Hmisc_3.17-1         
## [76] foreign_0.8-66        mgcv_1.8-10           survival_2.38-3      
## [79] RCurl_1.95-4.7        nnet_7.3-11           KernSmooth_2.23-15   
## [82] data.table_1.9.6      digest_0.6.9          diptest_0.75-7       
## [85] stats4_3.2.2          munsell_0.4.2         quadprog_1.5-5
```

