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
##  [1] earlywarnings_1.1.22 tseries_0.10-34      tgp_2.4-11          
##  [4] moments_0.14         gridExtra_2.0.0      RSQLite_1.0.0       
##  [7] DBI_0.3.1            googleVis_0.5.10     rdryad_0.2.0        
## [10] knitcitations_1.0.7  knitr_1.12           devtools_1.9.1      
## [13] limma_3.26.5         sorvi_0.7.35         ggplot2_2.0.0       
## [16] tidyr_0.3.1          dplyr_0.4.3          MASS_7.3-45         
## [19] netresponse_1.21.14  reshape2_1.4.1       mclust_5.1          
## [22] minet_3.28.0         Rgraphviz_2.14.0     graph_1.48.0        
## [25] phyloseq_1.14.0      microbiome_0.99.73   RPA_1.26.0          
## [28] affy_1.48.0          Biobase_2.30.0       BiocGenerics_0.16.1 
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      rjson_0.2.15          dynamicTreeCut_1.62  
##  [4] som_0.3-5             qvalue_2.2.2          XVector_0.10.0       
##  [7] affyio_1.40.0         AnnotationDbi_1.32.3  mvtnorm_1.0-3        
## [10] lubridate_1.5.0       RefManageR_0.10.5     xml2_0.1.2           
## [13] codetools_0.2-14      splines_3.2.2         doParallel_1.0.10    
## [16] impute_1.44.0         ade4_1.7-3            spam_1.3-0           
## [19] Formula_1.2-1         WGCNA_1.49            cluster_2.0.3        
## [22] GO.db_3.2.2           Kendall_2.2           oai_0.1.0            
## [25] httr_1.0.0            assertthat_0.1        Matrix_1.2-3         
## [28] lazyeval_0.1.10       formatR_1.2.1         acepack_1.3-3.3      
## [31] tools_3.2.2           igraph_1.0.1          gtable_0.1.2         
## [34] maps_3.0.2            Rcpp_0.12.3           Biostrings_2.38.3    
## [37] RJSONIO_1.3-0         multtest_2.26.0       biom_0.3.12          
## [40] ape_3.4               preprocessCore_1.32.0 nlme_3.1-122         
## [43] iterators_1.0.8       lmtest_0.9-34         fastcluster_1.1.16   
## [46] stringr_1.0.0         XML_3.98-1.3          zlibbioc_1.16.0      
## [49] zoo_1.7-12            scales_0.3.0          BiocInstaller_1.20.1 
## [52] solr_0.1.6            RColorBrewer_1.1-2    fields_8.3-6         
## [55] curl_0.9.5            memoise_0.2.1         rpart_4.1-10         
## [58] latticeExtra_0.6-26   stringi_1.0-1         maptree_1.4-7        
## [61] highr_0.5.1           S4Vectors_0.8.7       foreach_1.4.3        
## [64] nortest_1.0-4         permute_0.8-4         boot_1.3-17          
## [67] bibtex_0.4.0          chron_2.3-47          bitops_1.0-6         
## [70] matrixStats_0.50.1    dmt_0.8.20            evaluate_0.8         
## [73] lattice_0.20-33       labeling_0.3          plyr_1.8.3           
## [76] magrittr_1.5          R6_2.1.2              IRanges_2.4.6        
## [79] Hmisc_3.17-1          foreign_0.8-66        mgcv_1.8-10          
## [82] survival_2.38-3       RCurl_1.95-4.7        nnet_7.3-11          
## [85] KernSmooth_2.23-15    data.table_1.9.6      vegan_2.3-3          
## [88] digest_0.6.9          diptest_0.75-7        stats4_3.2.2         
## [91] munsell_0.4.2         quadprog_1.5-5
```

