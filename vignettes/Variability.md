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
##  [1] intergraph_2.0-2    sna_2.3-2           network_1.13.0     
##  [4] ggnet_0.1.0         GGally_1.0.1        knitcitations_1.0.7
##  [7] knitr_1.12          devtools_1.9.1      limma_3.26.5       
## [10] sorvi_0.7.35        ggplot2_2.0.0       tidyr_0.3.1        
## [13] dplyr_0.4.3         MASS_7.3-45         netresponse_1.21.14
## [16] reshape2_1.4.1      mclust_5.1          minet_3.28.0       
## [19] Rgraphviz_2.14.0    graph_1.48.0        phyloseq_1.14.0    
## [22] microbiome_0.99.73  RPA_1.26.0          affy_1.48.0        
## [25] Biobase_2.30.0      BiocGenerics_0.16.1
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-122          bitops_1.0-6          solr_0.1.6           
##  [4] lubridate_1.5.0       oai_0.1.0             RColorBrewer_1.1-2   
##  [7] httr_1.0.0            tools_3.2.2           R6_2.1.2             
## [10] vegan_2.3-3           affyio_1.40.0         rpart_4.1-10         
## [13] KernSmooth_2.23-15    dmt_0.8.20            lazyeval_0.1.10      
## [16] nortest_1.0-4         DBI_0.3.1             mgcv_1.8-10          
## [19] colorspace_1.2-6      permute_0.8-4         ade4_1.7-3           
## [22] moments_0.14          curl_0.9.5            preprocessCore_1.32.0
## [25] chron_2.3-47          rdryad_0.2.0          formatR_1.2.1        
## [28] xml2_0.1.2            labeling_0.3          tseries_0.10-34      
## [31] diptest_0.75-7        scales_0.3.0          lmtest_0.9-34        
## [34] mvtnorm_1.0-3         quadprog_1.5-5        tgp_2.4-11           
## [37] stringr_1.0.0         digest_0.6.9          earlywarnings_1.1.22 
## [40] XVector_0.10.0        bibtex_0.4.0          highr_0.5.1          
## [43] maps_3.0.2            BiocInstaller_1.20.1  zoo_1.7-12           
## [46] RCurl_1.95-4.7        magrittr_1.5          Matrix_1.2-3         
## [49] Rcpp_0.12.3           munsell_0.4.2         S4Vectors_0.8.7      
## [52] maptree_1.4-7         ape_3.4               RefManageR_0.10.5    
## [55] stringi_1.0-1         RJSONIO_1.3-0         zlibbioc_1.16.0      
## [58] plyr_1.8.3            qvalue_2.2.2          lattice_0.20-33      
## [61] Biostrings_2.38.3     splines_3.2.2         multtest_2.26.0      
## [64] igraph_1.0.1          boot_1.3-17           rjson_0.2.15         
## [67] codetools_0.2-14      stats4_3.2.2          XML_3.98-1.3         
## [70] evaluate_0.8          biom_0.3.12           data.table_1.9.6     
## [73] spam_1.3-0            foreach_1.4.3         gtable_0.1.2         
## [76] reshape_0.8.5         assertthat_0.1        Kendall_2.2          
## [79] survival_2.38-3       iterators_1.0.8       som_0.3-5            
## [82] memoise_0.2.1         IRanges_2.4.6         fields_8.3-6         
## [85] cluster_2.0.3
```

