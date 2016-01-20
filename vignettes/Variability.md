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

```
## Error in eval(expr, envir, enclos): could not find function "sample_data"
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

```
## Error in estimate_variability(x, "intraindividual"): The following variables needed by estimate_variability function type=intraindividual are missing from sample metadata: time
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

![plot of chunk variability-intra](figure/variability-intra-1.png)


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
##  [1] RSQLite_1.0.0       DBI_0.3.1           microbiome_0.99.71 
##  [4] devtools_1.9.1      sorvi_0.7.34        ggplot2_2.0.0      
##  [7] netresponse_1.21.14 reshape2_1.4.1      mclust_5.1         
## [10] minet_3.28.0        Rgraphviz_2.14.0    graph_1.48.0       
## [13] dplyr_0.4.3         limma_3.26.5        googleVis_0.5.10   
## [16] RPA_1.26.0          affy_1.48.0         Biobase_2.30.0     
## [19] BiocGenerics_0.16.1 knitcitations_1.0.7 knitr_1.12         
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.2-6      rjson_0.2.15          dynamicTreeCut_1.62  
##   [4] som_0.3-5             qvalue_2.2.2          corpcor_1.6.8        
##   [7] XVector_0.10.0        affyio_1.40.0         AnnotationDbi_1.32.3 
##  [10] mvtnorm_1.0-3         lubridate_1.5.0       RefManageR_0.10.5    
##  [13] xml2_0.1.2            codetools_0.2-14      splines_3.2.2        
##  [16] doParallel_1.0.10     impute_1.44.0         mixOmics_5.2.0       
##  [19] tgp_2.4-11            ade4_1.7-3            spam_1.3-0           
##  [22] Formula_1.2-1         phyloseq_1.14.0       WGCNA_1.49           
##  [25] cluster_2.0.3         GO.db_3.2.2           Kendall_2.2          
##  [28] oai_0.1.0             httr_1.0.0            lazyeval_0.1.10      
##  [31] assertthat_0.1        Matrix_1.2-3          formatR_1.2.1        
##  [34] acepack_1.3-3.3       tools_3.2.2           igraph_1.0.1         
##  [37] rdryad_0.2.0          gtable_0.1.2          maps_3.0.2           
##  [40] Rcpp_0.12.3           Biostrings_2.38.3     RJSONIO_1.3-0        
##  [43] multtest_2.26.0       biom_0.3.12           ape_3.4              
##  [46] preprocessCore_1.32.0 nlme_3.1-122          iterators_1.0.8      
##  [49] lmtest_0.9-34         fastcluster_1.1.16    stringr_1.0.0        
##  [52] XML_3.98-1.3          zlibbioc_1.16.0       MASS_7.3-45          
##  [55] zoo_1.7-12            scales_0.3.0          BiocInstaller_1.20.1 
##  [58] solr_0.1.6            RColorBrewer_1.1-2    fields_8.3-6         
##  [61] curl_0.9.4            gridExtra_2.0.0       memoise_0.2.1        
##  [64] rpart_4.1-10          latticeExtra_0.6-26   stringi_1.0-1        
##  [67] maptree_1.4-7         highr_0.5.1           S4Vectors_0.8.7      
##  [70] tseries_0.10-34       foreach_1.4.3         nortest_1.0-4        
##  [73] permute_0.8-4         boot_1.3-17           bibtex_0.4.0         
##  [76] chron_2.3-47          moments_0.14          bitops_1.0-6         
##  [79] matrixStats_0.50.1    rgl_0.95.1441         dmt_0.8.20           
##  [82] evaluate_0.8          lattice_0.20-33       labeling_0.3         
##  [85] plyr_1.8.3            magrittr_1.5          R6_2.1.1             
##  [88] IRanges_2.4.6         earlywarnings_1.1.22  Hmisc_3.17-1         
##  [91] foreign_0.8-66        mgcv_1.8-10           nnet_7.3-11          
##  [94] survival_2.38-3       RCurl_1.95-4.7        KernSmooth_2.23-15   
##  [97] ellipse_0.3-8         data.table_1.9.6      vegan_2.3-3          
## [100] digest_0.6.9          diptest_0.75-7        stats4_3.2.2         
## [103] munsell_0.4.2         quadprog_1.5-5
```

