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

Sarle's bimodality coefficient (see help(coefficient_of_bimodality)):


```r
bimodality.sarle <- bimodality(pseq0, method = "Sarle.finite.sample")
```


DIP test for multimodality (from the [diptest](https://cran.r-project.org/web/packages/diptest/index.html) package):


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
library(ggplot2)
grid.arrange(p1, p2, nrow = 1)
```

![plot of chunk stability2](figure/stability2-1.png) 


## Comparing bimodality and intermediate stability

The analysis suggests that bimodal population distribution across individuals is often associated with instable intermediate abundances within individuals. The specific bi-stable groups in the upper left corner were suggested to constitute bistable tipping elements of the human intestinal microbiota in [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html):


```r
taxa <- taxa_names(pseq0)
df <- data.frame(group = taxa,
                 intermediate.stability = intermediate.stability[taxa],
		 bimodality = bimodality.pb[taxa])
theme_set(theme_bw(20))
p <- ggplot(df, aes(x = intermediate.stability, y = bimodality, label = group))
p <- p + geom_text(size = 3)
p
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png) 


### Version information


```r
sessionInfo()
```

```
## R version 3.2.1 (2015-06-18)
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
##  [1] ggplot2_1.0.1       gridExtra_0.9.1     microbiome_0.99.61 
##  [4] RPA_1.24.0          affy_1.46.1         Biobase_2.28.0     
##  [7] BiocGenerics_0.14.0 phyloseq_1.13.2     knitr_1.10.5       
## [10] scimapClient_0.2.1 
## 
## loaded via a namespace (and not attached):
##   [1] nlme_3.1-120              bitops_1.0-6             
##   [3] RColorBrewer_1.1-2        GenomeInfoDb_1.4.1       
##   [5] tools_3.2.1               R6_2.1.0                 
##   [7] vegan_2.3-0               affyio_1.36.0            
##   [9] rpart_4.1-9               KernSmooth_2.23-14       
##  [11] Hmisc_3.16-0              nortest_1.0-3            
##  [13] DBI_0.3.1                 mgcv_1.8-6               
##  [15] colorspace_1.2-6          permute_0.8-4            
##  [17] ade4_1.7-2                nnet_7.3-9               
##  [19] DESeq2_1.8.1              moments_0.14             
##  [21] preprocessCore_1.30.0     chron_2.3-47             
##  [23] rdryad_0.1.1              formatR_1.2              
##  [25] RGCCA_2.0                 labeling_0.3             
##  [27] tseries_0.10-34           diptest_0.75-7           
##  [29] scales_0.2.5              lmtest_0.9-34            
##  [31] genefilter_1.50.0         quadprog_1.5-5           
##  [33] OAIHarvester_0.1-7        tgp_2.4-11               
##  [35] stringr_1.0.0             digest_0.6.8             
##  [37] foreign_0.8-63            earlywarnings_1.1.22     
##  [39] XVector_0.8.0             maps_2.3-9               
##  [41] RSQLite_1.0.0             BiocInstaller_1.18.3     
##  [43] zoo_1.7-12                gtools_3.5.0             
##  [45] BiocParallel_1.2.9        dplyr_0.4.2              
##  [47] acepack_1.3-3.3           RCurl_1.95-4.6           
##  [49] magrittr_1.5              Formula_1.2-1            
##  [51] futile.logger_1.4.1       Matrix_1.2-1             
##  [53] Rcpp_0.11.6               munsell_0.4.2            
##  [55] S4Vectors_0.6.2           maptree_1.4-7            
##  [57] ape_3.3                   proto_0.3-10             
##  [59] stringi_0.5-5             MASS_7.3-41              
##  [61] RJSONIO_1.3-0             zlibbioc_1.14.0          
##  [63] plyr_1.8.3                gdata_2.16.1             
##  [65] lattice_0.20-31           Biostrings_2.36.1        
##  [67] splines_3.2.1             multtest_2.24.0          
##  [69] annotate_1.46.1           locfit_1.5-9.1           
##  [71] igraph_0.7.1              GenomicRanges_1.20.5     
##  [73] boot_1.3-16               mixOmics_5.0-4           
##  [75] geneplotter_1.46.0        reshape2_1.4.1           
##  [77] codetools_0.2-11          stats4_3.2.1             
##  [79] futile.options_1.0.0      XML_3.98-1.3             
##  [81] evaluate_0.7              RcppArmadillo_0.5.200.1.0
##  [83] latticeExtra_0.6-26       biom_0.3.12              
##  [85] lambda.r_1.1.7            data.table_1.9.4         
##  [87] spam_1.0-1                foreach_1.4.2            
##  [89] gtable_0.1.2              assertthat_0.1           
##  [91] xtable_1.7-4              sorvi_0.7.30             
##  [93] Kendall_2.2               survival_2.38-2          
##  [95] pheatmap_1.0.2            iterators_1.0.7          
##  [97] som_0.3-5                 AnnotationDbi_1.30.1     
##  [99] IRanges_2.2.5             fields_8.2-1             
## [101] cluster_2.0.2             rgl_0.95.1247
```

