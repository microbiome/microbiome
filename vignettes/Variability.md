### Inter-individual heterogeneity (heterogeneity within group of samples)

Assess 'inter-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between samples and their mean for a given samples vs phylotypes matrix. For the illustration, calculate inter-individual stability (heterogeneity) separately for Placebo and LGG groups.

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
## Error in eval(expr, envir, enclos): object 'correlation' not found
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
## Error in eval(expr, envir, enclos): object 'correlation' not found
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
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] ggplot2_2.0.0       pROC_1.8            dplyr_0.4.3        
##  [4] limma_3.26.5        googleVis_0.5.10    knitcitations_1.0.7
##  [7] knitr_1.12          phyloseq_1.14.0     microbiome_0.99.73 
## [10] RPA_1.26.0          affy_1.48.0         Biobase_2.30.0     
## [13] BiocGenerics_0.16.1
## 
## loaded via a namespace (and not attached):
##  [1] netresponse_1.21.14   nlme_3.1-122          bitops_1.0-6         
##  [4] solr_0.1.6            lubridate_1.5.0       oai_0.1.0            
##  [7] RColorBrewer_1.1-2    httr_1.0.0            Rgraphviz_2.14.0     
## [10] tools_3.2.2           R6_2.1.2              vegan_2.3-3          
## [13] affyio_1.40.0         rpart_4.1-10          KernSmooth_2.23-15   
## [16] dmt_0.8.20            lazyeval_0.1.10       nortest_1.0-4        
## [19] DBI_0.3.1             mgcv_1.8-10           colorspace_1.2-6     
## [22] permute_0.8-4         ade4_1.7-3            moments_0.14         
## [25] preprocessCore_1.32.0 chron_2.3-47          rdryad_0.2.0         
## [28] graph_1.48.0          formatR_1.2.1         xml2_0.1.2           
## [31] labeling_0.3          tseries_0.10-34       diptest_0.75-7       
## [34] scales_0.3.0          lmtest_0.9-34         mvtnorm_1.0-3        
## [37] quadprog_1.5-5        tgp_2.4-11            stringr_1.0.0        
## [40] digest_0.6.9          earlywarnings_1.1.22  XVector_0.10.0       
## [43] bibtex_0.4.0          maps_3.0.2            BiocInstaller_1.20.1 
## [46] zoo_1.7-12            mclust_5.1            RCurl_1.95-4.7       
## [49] magrittr_1.5          Matrix_1.2-3          Rcpp_0.12.3          
## [52] munsell_0.4.2         S4Vectors_0.8.7       maptree_1.4-7        
## [55] ape_3.4               RefManageR_0.10.5     minet_3.28.0         
## [58] stringi_1.0-1         MASS_7.3-45           RJSONIO_1.3-0        
## [61] zlibbioc_1.16.0       plyr_1.8.3            qvalue_2.2.2         
## [64] grid_3.2.2            lattice_0.20-33       Biostrings_2.38.3    
## [67] splines_3.2.2         multtest_2.26.0       igraph_1.0.1         
## [70] boot_1.3-17           rjson_0.2.15          reshape2_1.4.1       
## [73] codetools_0.2-14      stats4_3.2.2          XML_3.98-1.3         
## [76] evaluate_0.8          biom_0.3.12           data.table_1.9.6     
## [79] spam_1.3-0            foreach_1.4.3         gtable_0.1.2         
## [82] tidyr_0.3.1           assertthat_0.1        sorvi_0.7.35         
## [85] Kendall_2.2           survival_2.38-3       iterators_1.0.8      
## [88] som_0.3-5             IRanges_2.4.6         fields_8.3-6         
## [91] cluster_2.0.3
```

