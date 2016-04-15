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
##  [1] earlywarnings_1.1.22 tseries_0.10-34      tgp_2.4-14          
##  [4] moments_0.14         gridExtra_2.2.1      knitcitations_1.0.7 
##  [7] knitr_1.12.3         intergraph_2.0-2     sna_2.3-2           
## [10] network_1.13.0       ggnet_0.1.0          GGally_1.0.1        
## [13] devtools_1.11.0      limma_3.26.9         sorvi_0.7.40        
## [16] ggplot2_2.1.0        tidyr_0.4.1          dplyr_0.4.3         
## [19] MASS_7.3-45          netresponse_1.20.15  reshape2_1.4.1      
## [22] mclust_5.2           minet_3.28.0         Rgraphviz_2.14.0    
## [25] graph_1.48.0         microbiome_0.99.78   RPA_1.27.3          
## [28] phyloseq_1.14.0      affy_1.48.0          Biobase_2.30.0      
## [31] BiocGenerics_0.16.1 
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-126          bitops_1.0-6          lubridate_1.5.6      
##  [4] httr_1.1.0            RColorBrewer_1.1-2    tools_3.2.2          
##  [7] R6_2.1.2              vegan_2.3-5           affyio_1.40.0        
## [10] rpart_4.1-10          KernSmooth_2.23-15    dmt_0.8.20           
## [13] lazyeval_0.1.10       nortest_1.0-4         DBI_0.3.1            
## [16] mgcv_1.8-12           colorspace_1.2-6      permute_0.9-0        
## [19] ade4_1.7-4            withr_1.0.1           preprocessCore_1.32.0
## [22] chron_2.3-47          formatR_1.3           labeling_0.3         
## [25] diptest_0.75-7        scales_0.4.0          lmtest_0.9-34        
## [28] mvtnorm_1.0-5         quadprog_1.5-5        stringr_1.0.0        
## [31] digest_0.6.9          XVector_0.10.0        bibtex_0.4.0         
## [34] maps_3.1.0            BiocInstaller_1.20.1  zoo_1.7-12           
## [37] RCurl_1.95-4.8        magrittr_1.5          Matrix_1.2-4         
## [40] Rcpp_0.12.4           munsell_0.4.3         S4Vectors_0.8.11     
## [43] maptree_1.4-7         ape_3.4               RefManageR_0.10.13   
## [46] stringi_1.0-1         RJSONIO_1.3-0         zlibbioc_1.16.0      
## [49] plyr_1.8.3            qvalue_2.2.2          lattice_0.20-33      
## [52] Biostrings_2.38.4     splines_3.2.2         multtest_2.26.0      
## [55] igraph_1.0.1          boot_1.3-18           codetools_0.2-14     
## [58] stats4_3.2.2          XML_3.98-1.4          evaluate_0.8.3       
## [61] biom_0.3.12           data.table_1.9.6      spam_1.3-0           
## [64] foreach_1.4.3         gtable_0.2.0          reshape_0.8.5        
## [67] assertthat_0.1        Kendall_2.2           survival_2.38-3      
## [70] iterators_1.0.8       som_0.3-5             memoise_1.0.0        
## [73] IRanges_2.4.8         fields_8.3-6          cluster_2.0.3
```

