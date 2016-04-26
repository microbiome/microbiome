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

### Version information


```r
sessionInfo()
```

```
## R version 3.2.3 (2015-12-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04 LTS
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
##  [1] splines   stats4    grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] microbiome_0.99.81   vegan_2.3-5          permute_0.9-0       
##  [4] earlywarnings_1.1.22 tseries_0.10-34      tgp_2.4-14          
##  [7] moments_0.14         gridExtra_2.2.1      scales_0.4.0        
## [10] SpiecEasi_0.1        VGAM_1.0-0           huge_1.2.7          
## [13] igraph_1.0.1         lattice_0.20-33      Matrix_1.2-4        
## [16] knitcitations_1.0.7  knitr_1.12.3         intergraph_2.0-2    
## [19] sna_2.3-2            network_1.13.0       ggnet_0.1.0         
## [22] GGally_1.0.1         devtools_1.10.0      limma_3.26.9        
## [25] sorvi_0.7.43         tibble_1.0           ggplot2_2.1.0       
## [28] tidyr_0.4.1          dplyr_0.4.3          MASS_7.3-45         
## [31] netresponse_1.20.15  reshape2_1.4.1       mclust_5.2          
## [34] minet_3.28.0         Rgraphviz_2.14.0     graph_1.48.0        
## [37] phyloseq_1.14.0     
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-126        bitops_1.0-6        lubridate_1.5.0    
##  [4] RColorBrewer_1.1-2  httr_1.1.0          tools_3.2.3        
##  [7] R6_2.1.2            rpart_4.1-10        KernSmooth_2.23-15 
## [10] dmt_0.8.20          lazyeval_0.1.10     nortest_1.0-4      
## [13] DBI_0.3.1           BiocGenerics_0.16.1 mgcv_1.8-12        
## [16] colorspace_1.2-6    ade4_1.7-4          withr_1.0.1        
## [19] chron_2.3-47        Biobase_2.30.0      formatR_1.3        
## [22] labeling_0.3        lmtest_0.9-34       mvtnorm_1.0-5      
## [25] quadprog_1.5-5      stringr_1.0.0       digest_0.6.9       
## [28] XVector_0.10.0      bibtex_0.4.0        highr_0.5.1        
## [31] maps_3.1.0          zoo_1.7-12          RCurl_1.95-4.8     
## [34] magrittr_1.5        Rcpp_0.12.4         munsell_0.4.3      
## [37] S4Vectors_0.8.11    maptree_1.4-7       ape_3.4            
## [40] RefManageR_0.10.6   stringi_1.0-1       RJSONIO_1.3-0      
## [43] zlibbioc_1.16.0     plyr_1.8.3          qvalue_2.2.2       
## [46] parallel_3.2.3      Biostrings_2.38.4   multtest_2.26.0    
## [49] boot_1.3-18         codetools_0.2-14    XML_3.98-1.4       
## [52] evaluate_0.8.3      biom_0.3.12         data.table_1.9.6   
## [55] spam_1.3-0          foreach_1.4.3       gtable_0.2.0       
## [58] reshape_0.8.5       assertthat_0.1      roxygen2_5.0.1     
## [61] Kendall_2.2         survival_2.38-3     iterators_1.0.8    
## [64] som_0.3-5           memoise_1.0.0       IRanges_2.4.8      
## [67] fields_8.3-6        cluster_2.0.3
```

