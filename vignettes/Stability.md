# Microbiome stability analysis

Download example data - HITChip Atlas of 130 genus-like taxa across 1006 healthy western adults from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). A subset of 76 subjects have also short time series available for temporal stability analysis:


```r
library(microbiome)
pseq <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```r
# Let us keep only prevalent taxa
# (HITChip signal >3 in >20 percent of the samples)
pseq <- filter_prevalent(pseq, detection.threshold = 10^3,
     			       prevalence.threshold = 0.2)
```



## Quantify intermediate stability 

It has been reported that certain microbial groups exhibit bi-stable
abundance distributions with distinct peaks at low and high
abundances, and an instable intermediate abundance range.

Instability at the intermediate abundance range is hence one indicator
of bi-stability.

[Lahti et
al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html))
use straightforward correlation analysis to quantify how the distance
from the intermediate abundance region (50% quantile) is associated
with the observed shifts between consecutive time points. This can be
calculated with:


```r
intermediate.stability <- intermediate_stability(pseq, output = "scores")
```


## Quantify bimodality 

Bistable systems often exhibit bimodal population patterns. Hence,
bimodality of the abundance distribution provides another indicator of
bistability, although other explanations for bimodality (sampling
biases; heterogeneous population structure etc.) are also possible.

Multiple bimodality scores are available. Perform cross-sectional
analysis of bimodality by including only the samples from the zero
time point:


```r
pseq0 <- subset_samples(pseq, time == 0 & DNA_extraction_method == "r")
```


[Multimodality score with potential analysis with
bootstrap](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)



```r
bimodality.pb <- bimodality(pseq0, method = "potential_bootstrap")
```

Sarle's bimodality coefficient (see help(coefficient_of_bimodality)):


```r
bimodality.sarle <- bimodality(pseq0, method = "Sarle.finite.sample")
```


DIP test for multimodality (not included in the microbiome package):


```r
bimodality.dip <- bimodality(pseq0, method = "dip")
```


Visualize population densities 


```r
# Pick the most and least bimodal taxa as examples
bimodality <- bimodality.pb
unimodal <- names(which.min(bimodality))
bimodal  <- names(which.max(bimodality))

# Visualize population frequencies
p1 <- plot_density(pseq, otu.name = unimodal, log10 = TRUE) 
p2 <- plot_density(pseq, otu.name = bimodal,  log10 = TRUE) 
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
```

```
## Error: could not find function "ggplotGrob"
```


## Comparing bimodality and intermediate stability

The analysis suggests that bimodal population distribution across individuals is often associated with instable intermediate abundances within individuals. The specific bi-stable groups in the upper left corner were suggested in [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html):


```r
df <- data.frame(group = tax,
                 intermediate.stability = intermediate.stability[tax],
		 bimodality = bimodality.pb[tax])
```

```
## Error in data.frame(group = tax, intermediate.stability = intermediate.stability[tax], : object 'tax' not found
```

```r
theme_set(theme_bw(20))
```

```
## Error in eval(expr, envir, enclos): could not find function "theme_set"
```

```r
p <- ggplot(df, aes(x = intermediate.stability, y = bimodality, label = group))
```

```
## Error in eval(expr, envir, enclos): could not find function "ggplot"
```

```r
p <- p + geom_text(size = 3)
```

```
## Error in eval(expr, envir, enclos): could not find function "geom_text"
```

```r
p
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png) 


### TODO

 * tipping element heatmap

 * stability line plot

 * classification of the taxa into the distinct abundance types

### Version information


```r
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-unknown-linux-gnu (64-bit)
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
## [1] gridExtra_0.9.1     microbiome_0.99.52  RPA_1.24.0         
## [4] affy_1.46.0         Biobase_2.28.0      BiocGenerics_0.14.0
## [7] phyloseq_1.13.2     knitr_1.10.5        scimapClient_0.2.1 
## 
## loaded via a namespace (and not attached):
##   [1] nlme_3.1-120              bitops_1.0-6             
##   [3] RColorBrewer_1.1-2        GenomeInfoDb_1.4.0       
##   [5] tools_3.2.0               vegan_2.2-1              
##   [7] affyio_1.36.0             rpart_4.1-9              
##   [9] KernSmooth_2.23-14        Hmisc_3.16-0             
##  [11] nortest_1.0-3             DBI_0.3.1                
##  [13] mgcv_1.8-6                colorspace_1.2-6         
##  [15] permute_0.8-4             ade4_1.7-2               
##  [17] nnet_7.3-9                DESeq2_1.8.1             
##  [19] moments_0.14              preprocessCore_1.30.0    
##  [21] chron_2.3-45              rdryad_0.1.1             
##  [23] formatR_1.2               RGCCA_2.0                
##  [25] tseries_0.10-34           diptest_0.75-6           
##  [27] scales_0.2.4              lmtest_0.9-33            
##  [29] genefilter_1.50.0         quadprog_1.5-5           
##  [31] OAIHarvester_0.1-7        tgp_2.4-11               
##  [33] stringr_1.0.0             digest_0.6.8             
##  [35] foreign_0.8-63            earlywarnings_1.1.19     
##  [37] XVector_0.8.0             maps_2.3-9               
##  [39] RSQLite_1.0.0             BiocInstaller_1.18.2     
##  [41] zoo_1.7-12                gtools_3.4.2             
##  [43] BiocParallel_1.2.1        dplyr_0.4.1              
##  [45] acepack_1.3-3.3           RCurl_1.95-4.6           
##  [47] magrittr_1.5              Formula_1.2-1            
##  [49] futile.logger_1.4.1       Matrix_1.2-0             
##  [51] Rcpp_0.11.6               munsell_0.4.2            
##  [53] S4Vectors_0.6.0           maptree_1.4-7            
##  [55] ape_3.2                   proto_0.3-10             
##  [57] stringi_0.4-1             MASS_7.3-40              
##  [59] RJSONIO_1.3-0             zlibbioc_1.14.0          
##  [61] plyr_1.8.2                gdata_2.16.1             
##  [63] lattice_0.20-31           Biostrings_2.36.1        
##  [65] splines_3.2.0             multtest_2.24.0          
##  [67] annotate_1.46.0           locfit_1.5-9.1           
##  [69] igraph_0.7.1              GenomicRanges_1.20.3     
##  [71] boot_1.3-16               mixOmics_5.0-4           
##  [73] geneplotter_1.46.0        reshape2_1.4.1           
##  [75] codetools_0.2-11          stats4_3.2.0             
##  [77] futile.options_1.0.0      XML_3.98-1.1             
##  [79] evaluate_0.7              RcppArmadillo_0.5.100.1.0
##  [81] latticeExtra_0.6-26       biom_0.3.12              
##  [83] lambda.r_1.1.7            data.table_1.9.4         
##  [85] spam_1.0-1                foreach_1.4.2            
##  [87] gtable_0.1.2              assertthat_0.1           
##  [89] ggplot2_1.0.1             xtable_1.7-4             
##  [91] sorvi_0.7.23              Kendall_2.2              
##  [93] survival_2.38-1           pheatmap_1.0.2           
##  [95] iterators_1.0.7           som_0.3-5                
##  [97] AnnotationDbi_1.30.1      IRanges_2.2.1            
##  [99] fields_8.2-1              cluster_2.0.1            
## [101] rgl_0.95.1247
```

