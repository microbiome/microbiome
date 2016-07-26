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


## Further visualization tools

Draw regression curve with smoothed error bars based on
the [Visually-Weighted Regression](http://www.fight-entropy.com/2012/07/visually-weighted-regression.html) by Solomon M. Hsiang. The sorvi implementation extends [Felix Schonbrodt's original code](http://www.nicebread.de/visually-weighted-watercolor-plots-new-variants-please-vote/).


```r
data(atlas1006)
p <- plot_regression(diversity ~ age, sample_data(atlas1006))
print(p)
```

![plot of chunk variability-regression](figure/variability-regression-1.png)

### Version information


```r
sessionInfo()
```

```
## R version 3.3.1 (2016-06-21)
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
##  [1] tcltk     grid      parallel  stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] earlywarnings_1.1.22  tseries_0.10-35       tgp_2.4-14           
##  [4] moments_0.14          ggrepel_0.5           gridExtra_2.2.1      
##  [7] HITChipDB_0.6.31      RPA_1.27.41           affy_1.50.0          
## [10] Biobase_2.32.0        RMySQL_0.10.9         preprocessCore_1.34.0
## [13] DBI_0.4-1             microbiome_0.99.86    FD_1.0-12            
## [16] geometry_0.3-6        magic_1.5-6           abind_1.4-3          
## [19] ape_3.5               ade4_1.7-4            RColorBrewer_1.1-2   
## [22] vegan_2.4-0           lattice_0.20-33       permute_0.9-0        
## [25] knitcitations_1.0.7   knitr_1.13            intergraph_2.0-2     
## [28] sna_2.3-2             network_1.13.0        ggnet_0.1.0          
## [31] GGally_1.2.0          devtools_1.12.0       limma_3.28.14        
## [34] sorvi_0.7.47          tibble_1.1            ggplot2_2.1.0        
## [37] tidyr_0.5.1           dplyr_0.5.0           MASS_7.3-45          
## [40] netresponse_1.3.17    reshape2_1.4.1        mclust_5.2           
## [43] minet_3.30.0          Rgraphviz_2.16.0      graph_1.50.0         
## [46] BiocGenerics_0.18.0   phyloseq_1.16.2      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-6      dynamicTreeCut_1.63-1 som_0.3-5            
##  [4] qvalue_2.4.2          XVector_0.12.0        affyio_1.42.0        
##  [7] AnnotationDbi_1.34.3  mvtnorm_1.0-5         lubridate_1.5.6      
## [10] RefManageR_0.10.13    codetools_0.2-14      splines_3.3.1        
## [13] doParallel_1.0.10     impute_1.46.0         spam_1.3-0           
## [16] Formula_1.2-1         jsonlite_1.0          Cairo_1.5-9          
## [19] WGCNA_1.51            cluster_2.0.4         GO.db_3.3.0          
## [22] Kendall_2.2           httr_1.2.1            assertthat_0.1       
## [25] Matrix_1.2-6          lazyeval_0.2.0        formatR_1.4          
## [28] acepack_1.3-3.3       tools_3.3.1           igraph_1.0.1         
## [31] gtable_0.2.0          maps_3.1.0            Rcpp_0.12.5          
## [34] Biostrings_2.40.2     RJSONIO_1.3-0         multtest_2.28.0      
## [37] nlme_3.1-128          iterators_1.0.8       lmtest_0.9-34        
## [40] fastcluster_1.1.20    stringr_1.0.0         XML_3.98-1.4         
## [43] zoo_1.7-13            zlibbioc_1.18.0       scales_0.4.0         
## [46] BiocInstaller_1.22.3  biomformat_1.0.2      rhdf5_2.16.0         
## [49] fields_8.4-1          curl_0.9.7            memoise_1.0.0        
## [52] rpart_4.1-10          reshape_0.8.5         latticeExtra_0.6-28  
## [55] stringi_1.1.1         maptree_1.4-7         RSQLite_1.0.0        
## [58] highr_0.6             S4Vectors_0.10.1      nortest_1.0-4        
## [61] foreach_1.4.3         boot_1.3-18           bibtex_0.4.0         
## [64] chron_2.3-47          matrixStats_0.50.2    bitops_1.0-6         
## [67] dmt_0.8.20            evaluate_0.9          labeling_0.3         
## [70] plyr_1.8.4            magrittr_1.5          R6_2.1.2             
## [73] IRanges_2.6.1         Hmisc_3.17-4          foreign_0.8-66       
## [76] withr_1.0.2           mgcv_1.8-12           survival_2.39-5      
## [79] RCurl_1.95-4.8        nnet_7.3-12           KernSmooth_2.23-15   
## [82] data.table_1.9.6      git2r_0.15.0          digest_0.6.9         
## [85] stats4_3.3.1          munsell_0.4.3         quadprog_1.5-5
```

