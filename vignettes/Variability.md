### Inter-individual variability (variability within group of samples)

Assess 'inter-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between samples and their mean for a given samples vs phylotypes matrix. For the illustration, calculate inter-individual stability (heterogeneity) separately for Placebo and LGG groups.

Load example data


```r
library(microbiome)
x <- download_microbiome("dietswap")
```

```
## Error in curl::curl_fetch_memory(url, handle = handle): Timeout was reached
```

```r
# Add time field (two time points needed within each group for the 
# intraindividual method)
sample_data(x)$time <- sample_data(x)$timepoint.within.group
```

```
## Error in access(object, "sam_data", errorIfNULL): sam_data slot is empty.
```


Heterogeneity across subjects within a group


```r
res <- estimate_heterogeneity(x, "interindividual")
```

```
## Error in eval(expr, envir, enclos): could not find function "estimate_heterogeneity"
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
```

```
## Error: ggplot2 doesn't know how to deal with data of class matrix
```

```r
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Inter-individual heterogeneity (p=", round(res$p.value, 6), ")"))
```

```
## Error in round(res$p.value, 6): non-numeric argument to mathematical function
```

```r
p <- p + ylab("Correlation")
print(p)
```

```
## Error: Aesthetics must be either length 1 or the same as the data (20): x, y, label
```

![plot of chunk heterogeneity-example2d](figure/heterogeneity-example2d-1.png)


### Intra-individual stability

Heterogeneity within subjects over time (intra-individual stability). Assess 'intra-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between two time points within subjects, calculated separately within each group. For illustration, check intra-individual stability (heterogeneity) separately for Placebo and LGG groups.


```r
res <- estimate_heterogeneity(x, "intraindividual")
```

```
## Error in eval(expr, envir, enclos): could not find function "estimate_heterogeneity"
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
```

```
## Error: ggplot2 doesn't know how to deal with data of class matrix
```

```r
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Intra-individual heterogeneity (p=", round(res$p.value, 6), ")"))
```

```
## Error in round(res$p.value, 6): non-numeric argument to mathematical function
```

```r
p <- p + ylab("Correlation")
print(p)
```

```
## Error: Aesthetics must be either length 1 or the same as the data (20): x, y, label
```

![plot of chunk heterogeneity-intra](figure/heterogeneity-intra-1.png)


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
##  [1] gridExtra_2.0.0     rdryad_0.2.0        pROC_1.8           
##  [4] sorvi_0.7.35        RSQLite_1.0.0       DBI_0.3.1          
##  [7] ggplot2_2.0.0       dplyr_0.4.3         limma_3.26.5       
## [10] googleVis_0.5.10    netresponse_1.21.14 reshape2_1.4.1     
## [13] mclust_5.1          minet_3.28.0        Rgraphviz_2.14.0   
## [16] graph_1.48.0        knitcitations_1.0.7 knitr_1.12         
## [19] phyloseq_1.14.0     microbiome_0.99.73  RPA_1.26.0         
## [22] affy_1.48.0         Biobase_2.30.0      BiocGenerics_0.16.1
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      rjson_0.2.15          dynamicTreeCut_1.62  
##  [4] som_0.3-5             qvalue_2.2.2          XVector_0.10.0       
##  [7] affyio_1.40.0         AnnotationDbi_1.32.3  mvtnorm_1.0-3        
## [10] lubridate_1.5.0       RefManageR_0.10.5     xml2_0.1.2           
## [13] codetools_0.2-14      splines_3.2.2         doParallel_1.0.10    
## [16] impute_1.44.0         tgp_2.4-11            ade4_1.7-3           
## [19] spam_1.3-0            Formula_1.2-1         WGCNA_1.49           
## [22] cluster_2.0.3         GO.db_3.2.2           Kendall_2.2          
## [25] oai_0.1.0             httr_1.0.0            assertthat_0.1       
## [28] Matrix_1.2-3          lazyeval_0.1.10       formatR_1.2.1        
## [31] acepack_1.3-3.3       tools_3.2.2           igraph_1.0.1         
## [34] gtable_0.1.2          maps_3.0.2            Rcpp_0.12.3          
## [37] Biostrings_2.38.3     RJSONIO_1.3-0         multtest_2.26.0      
## [40] biom_0.3.12           ape_3.4               preprocessCore_1.32.0
## [43] nlme_3.1-122          iterators_1.0.8       lmtest_0.9-34        
## [46] fastcluster_1.1.16    stringr_1.0.0         XML_3.98-1.3         
## [49] zlibbioc_1.16.0       MASS_7.3-45           zoo_1.7-12           
## [52] scales_0.3.0          BiocInstaller_1.20.1  solr_0.1.6           
## [55] RColorBrewer_1.1-2    fields_8.3-6          curl_0.9.5           
## [58] rpart_4.1-10          latticeExtra_0.6-26   stringi_1.0-1        
## [61] maptree_1.4-7         highr_0.5.1           S4Vectors_0.8.7      
## [64] tseries_0.10-34       foreach_1.4.3         nortest_1.0-4        
## [67] permute_0.8-4         boot_1.3-17           bibtex_0.4.0         
## [70] chron_2.3-47          moments_0.14          bitops_1.0-6         
## [73] matrixStats_0.50.1    dmt_0.8.20            evaluate_0.8         
## [76] lattice_0.20-33       labeling_0.3          plyr_1.8.3           
## [79] magrittr_1.5          R6_2.1.2              IRanges_2.4.6        
## [82] earlywarnings_1.1.22  Hmisc_3.17-1          foreign_0.8-66       
## [85] mgcv_1.8-10           survival_2.38-3       RCurl_1.95-4.7       
## [88] nnet_7.3-11           KernSmooth_2.23-15    data.table_1.9.6     
## [91] vegan_2.3-3           digest_0.6.9          diptest_0.75-7       
## [94] tidyr_0.3.1           stats4_3.2.2          munsell_0.4.2        
## [97] quadprog_1.5-5
```

