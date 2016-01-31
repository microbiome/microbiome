# Microbiome stability analysis

Download example data - HITChip Atlas of 130 genus-like taxa across 1006 healthy western adults from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). A subset of 76 subjects have also short time series available for temporal stability analysis:


```r
library(microbiome)
pseq <- download_microbiome("atlas1006")
```

```
## Error in curl::curl_fetch_memory(url, handle = handle): Timeout was reached
```

```r
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

```
## Error in eval(expr, envir, enclos): object 'DNA_extraction_method' not found
```


Multimodality score using [potential analysis with
bootstrap](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)



```r
bimodality.pb <- bimodality(pseq0, method = "potential_bootstrap")
```

```
## Error in is.vector(x): object 'pseq0' not found
```

Sarle's bimodality coefficient (see help(coefficient_of_bimodality)):


```r
bimodality.sarle <- bimodality(pseq0, method = "Sarle.finite.sample")
```

```
## Error in is.vector(x): object 'pseq0' not found
```


DIP test for multimodality (from the [diptest](https://cran.r-project.org/web/packages/diptest/index.html) package):


```r
bimodality.dip <- bimodality(pseq0, method = "dip")
```

```
## Error in is.vector(x): object 'pseq0' not found
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
## Error in which.min(bimodality): cannot coerce type 'closure' to vector of type 'double'
```

```r
bimodal  <- names(which.max(bimodality))
```

```
## Error in which.max(bimodality): cannot coerce type 'closure' to vector of type 'double'
```

```r
# Visualize population frequencies
p1 <- plot_density(pseq, otu.name = unimodal, log10 = TRUE) 
```

```
## Error in plot_density(pseq, otu.name = unimodal, log10 = TRUE): unused argument (otu.name = unimodal)
```

```r
p2 <- plot_density(pseq, otu.name = bimodal,  log10 = TRUE) 
```

```
## Error in plot_density(pseq, otu.name = bimodal, log10 = TRUE): unused argument (otu.name = bimodal)
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
```

```
## Error in taxa_names(pseq0): error in evaluating the argument 'physeq' in selecting a method for function 'taxa_names': Error: object 'pseq0' not found
```

```r
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
## Error: Aesthetics must be either length 1 or the same as the data (20): x, y, label
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png)


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

