# Microbiome stability analysis

Download example data - HITChip Atlas of 130 genus-like taxa across 1006 healthy western adults from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). A subset of 76 subjects have also short time series available for temporal stability analysis:


```r
library(microbiome)
pseq <- download_microbiome("atlas1006")

# Let us keep only prevalent taxa
# (HITChip signal >3 in >20 percent of the samples)
pseq <- filter_prevalent(pseq, detection.threshold = 10^3, prevalence.threshold = 0.2)
```



## Quantify intermediate stability 

It has been reported that certain microbial groups exhibit bi-stable
abundance distributions with distinct peaks at low and high
abundances, and an instable intermediate abundance range. Instability
at the intermediate abundance range is hence one indicator of
bi-stability. [Lahti et
al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html))
use straightforward correlation analysis to quantify how the distance
from the intermediate abundance region (50% quantile) is associated
with the observed shifts between consecutive time points. This can be
calculated with:


```r
intermediate.stability <- intermediate_stability(pseq, output = "scores")
```


## Quantify bimodality 

Bimodality of the abundance distribution provides another (indirect)
indicator of bistability, although other explanations such as sampling
biases etc. should be controlled. Multiple bimodality scores are
available.

Let is stick to cross-sectional analysis of bimodality and include
only the samples from the zero time point:


```r
pseq0 <- subset_samples(pseq, time == 0 & DNA_extraction_method == "r")
```


Multimodality score using [potential analysis with
bootstrap](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)



```r
bimodality.pb <- bimodality(pseq0, method = "potential_bootstrap")
```

```
## Error in eval(expr, envir, enclos): could not find function "bimodality"
```

Sarle's bimodality coefficient (see help(coefficient_of_bimodality)):


```r
bimodality.sarle <- bimodality(pseq0, method = "Sarle.finite.sample")
```

```
## Error in eval(expr, envir, enclos): could not find function "bimodality"
```


DIP test for multimodality (from the [diptest](https://cran.r-project.org/web/packages/diptest/index.html) package):


```r
bimodality.dip <- bimodality(pseq0, method = "dip")
```

```
## Error in eval(expr, envir, enclos): could not find function "bimodality"
```


Visualize population densities 


```r
# Pick the most and least bimodal taxa as examples
bimodality <- bimodality.pb
```

```
## Error in eval(expr, envir, enclos): object 'bimodality.pb' not found
```

```r
unimodal <- names(which.min(bimodality))
```

```
## Error in which.min(bimodality): object 'bimodality' not found
```

```r
bimodal  <- names(which.max(bimodality))
```

```
## Error in which.max(bimodality): object 'bimodality' not found
```

```r
# Visualize population frequencies
p1 <- plot_density(pseq, otu.name = unimodal, log10 = TRUE) 
```

```
## Error in plot_density(pseq, otu.name = unimodal, log10 = TRUE): object 'unimodal' not found
```

```r
p2 <- plot_density(pseq, otu.name = bimodal,  log10 = TRUE) 
```

```
## Error in plot_density(pseq, otu.name = bimodal, log10 = TRUE): object 'bimodal' not found
```

```r
library(gridExtra)
library(ggplot2)
grid.arrange(p1, p2, nrow = 1)
```

```
## Error in arrangeGrob(...): object 'p1' not found
```


## Comparing bimodality and intermediate stability

The analysis suggests that bimodal population distribution across individuals is often associated with instable intermediate abundances within individuals. The specific bi-stable groups in the upper left corner were suggested to constitute bistable tipping elements of the human intestinal microbiota in [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html):


```r
taxa <- taxa_names(pseq0)
df <- data.frame(group = taxa,
                 intermediate.stability = intermediate.stability[taxa],
		 bimodality = bimodality.pb[taxa])
```

```
## Error in data.frame(group = taxa, intermediate.stability = intermediate.stability[taxa], : object 'bimodality.pb' not found
```

```r
theme_set(theme_bw(20))
p <- ggplot(df, aes(x = intermediate.stability, y = bimodality, label = group))
p <- p + geom_text(size = 3)
p
```

```
## Error in eval(expr, envir, enclos): object 'bimodality' not found
```


### Version information


```r
sessionInfo()
```

```
## R version 3.2.2 (2015-08-14)
## Platform: x86_64-pc-linux-gnu (64-bit)
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
##  [1] googleVis_0.5.10    gridExtra_2.0.0     limma_3.24.15      
##  [4] mgcv_1.8-8          nlme_3.1-122        dplyr_0.4.3        
##  [7] netresponse_1.18.0  reshape_0.8.5       mclust_5.1         
## [10] minet_3.26.0        Rgraphviz_2.12.0    graph_1.46.0       
## [13] ggplot2_1.0.1       sorvi_0.7.32        microbiome_0.99.63 
## [16] RPA_1.24.0          affy_1.46.1         Biobase_2.28.0     
## [19] BiocGenerics_0.14.0 phyloseq_1.12.2     rdryad_0.1.1       
## [22] knitcitations_1.0.7 knitr_1.11         
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.2-6          qvalue_2.0.0             
##   [3] som_0.3-5                 futile.logger_1.4.1      
##   [5] XVector_0.8.0             OAIHarvester_0.1-7       
##   [7] RcppArmadillo_0.6.100.0.0 GenomicRanges_1.20.8     
##   [9] affyio_1.36.0             mvtnorm_1.0-3            
##  [11] AnnotationDbi_1.30.1      lubridate_1.3.3          
##  [13] RefManageR_0.8.63         codetools_0.2-14         
##  [15] splines_3.2.2             geneplotter_1.46.0       
##  [17] mixOmics_5.1.2            tgp_2.4-11               
##  [19] ade4_1.7-2                spam_1.3-0               
##  [21] Formula_1.2-1             annotate_1.46.1          
##  [23] cluster_2.0.3             pheatmap_1.0.7           
##  [25] Kendall_2.2               httr_1.0.0               
##  [27] lazyeval_0.1.10           assertthat_0.1           
##  [29] Matrix_1.2-2              formatR_1.2.1            
##  [31] acepack_1.3-3.3           tools_3.2.2              
##  [33] igraph_1.0.1              gtable_0.1.2             
##  [35] reshape2_1.4.1            maps_3.0.0-2             
##  [37] Rcpp_0.12.1               Biostrings_2.36.4        
##  [39] RJSONIO_1.3-0             multtest_2.24.0          
##  [41] biom_0.3.12               gdata_2.17.0             
##  [43] ape_3.3                   preprocessCore_1.30.0    
##  [45] iterators_1.0.8           lmtest_0.9-34            
##  [47] stringr_1.0.0             proto_0.3-10             
##  [49] gtools_3.5.0              XML_3.98-1.3             
##  [51] zlibbioc_1.14.0           MASS_7.3-44              
##  [53] zoo_1.7-12                scales_0.3.0             
##  [55] BiocInstaller_1.18.5      lambda.r_1.1.7           
##  [57] RColorBrewer_1.1-2        fields_8.3-5             
##  [59] memoise_0.2.1             rpart_4.1-10             
##  [61] latticeExtra_0.6-26       stringi_1.0-1            
##  [63] maptree_1.4-7             RSQLite_1.0.0            
##  [65] highr_0.5.1               genefilter_1.50.0        
##  [67] S4Vectors_0.6.6           tseries_0.10-34          
##  [69] foreach_1.4.3             nortest_1.0-4            
##  [71] permute_0.8-4             boot_1.3-17              
##  [73] BiocParallel_1.2.22       bibtex_0.4.0             
##  [75] chron_2.3-47              GenomeInfoDb_1.4.3       
##  [77] moments_0.14              bitops_1.0-6             
##  [79] dmt_0.8.20                rgl_0.95.1367            
##  [81] evaluate_0.8              lattice_0.20-33          
##  [83] labeling_0.3              plyr_1.8.3               
##  [85] magrittr_1.5              DESeq2_1.8.2             
##  [87] R6_2.1.1                  IRanges_2.2.9            
##  [89] earlywarnings_1.1.22      Hmisc_3.17-0             
##  [91] DBI_0.3.1                 foreign_0.8-66           
##  [93] survival_2.38-3           RCurl_1.95-4.7           
##  [95] nnet_7.3-11               futile.options_1.0.0     
##  [97] KernSmooth_2.23-15        ellipse_0.3-8            
##  [99] locfit_1.5-9.1            data.table_1.9.6         
## [101] vegan_2.3-1               digest_0.6.8             
## [103] diptest_0.75-7            xtable_1.7-4             
## [105] stats4_3.2.2              munsell_0.4.2            
## [107] quadprog_1.5-5
```

