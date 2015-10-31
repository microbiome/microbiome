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
## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] gridExtra_2.0.0     sorvi_0.7.32        mgcv_1.8-8         
##  [4] nlme_3.1-122        RSQLite_1.0.0       DBI_0.3.1          
##  [7] dplyr_0.4.3         limma_3.24.15       ggplot2_1.0.1      
## [10] microbiome_0.99.65  RPA_1.24.0          affy_1.46.1        
## [13] Biobase_2.28.0      BiocGenerics_0.14.0 phyloseq_1.12.2    
## [16] knitcitations_1.0.7 knitr_1.11         
## 
## loaded via a namespace (and not attached):
##   [1] colorspace_1.2-6          dynamicTreeCut_1.62      
##   [3] mclust_5.1                qvalue_2.0.0             
##   [5] som_0.3-5                 futile.logger_1.4.1      
##   [7] XVector_0.8.0             OAIHarvester_0.1-7       
##   [9] RcppArmadillo_0.6.100.0.0 GenomicRanges_1.20.8     
##  [11] affyio_1.36.0             mvtnorm_1.0-3            
##  [13] AnnotationDbi_1.30.1      lubridate_1.3.3          
##  [15] RefManageR_0.8.63         codetools_0.2-14         
##  [17] splines_3.2.2             doParallel_1.0.10        
##  [19] impute_1.42.0             geneplotter_1.46.0       
##  [21] mixOmics_5.1.2            tgp_2.4-11               
##  [23] ade4_1.7-2                spam_1.3-0               
##  [25] Formula_1.2-1             WGCNA_1.47               
##  [27] annotate_1.46.1           GO.db_3.1.2              
##  [29] cluster_2.0.3             graph_1.46.0             
##  [31] pheatmap_1.0.7            Kendall_2.2              
##  [33] httr_1.0.0                lazyeval_0.1.10          
##  [35] assertthat_0.1            Matrix_1.2-2             
##  [37] formatR_1.2.1             acepack_1.3-3.3          
##  [39] tools_3.2.2               igraph_1.0.1             
##  [41] rdryad_0.1.1              gtable_0.1.2             
##  [43] reshape2_1.4.1            maps_3.0.0-2             
##  [45] Rcpp_0.12.1               Biostrings_2.36.4        
##  [47] RJSONIO_1.3-0             multtest_2.24.0          
##  [49] biom_0.3.12               gdata_2.17.0             
##  [51] ape_3.3                   preprocessCore_1.30.0    
##  [53] iterators_1.0.8           lmtest_0.9-34            
##  [55] fastcluster_1.1.16        stringr_1.0.0            
##  [57] proto_0.3-10              gtools_3.5.0             
##  [59] XML_3.98-1.3              zlibbioc_1.14.0          
##  [61] MASS_7.3-44               zoo_1.7-12               
##  [63] scales_0.3.0              minet_3.26.0             
##  [65] BiocInstaller_1.18.5      lambda.r_1.1.7           
##  [67] RColorBrewer_1.1-2        fields_8.3-5             
##  [69] memoise_0.2.1             rpart_4.1-10             
##  [71] reshape_0.8.5             latticeExtra_0.6-26      
##  [73] stringi_1.0-1             maptree_1.4-7            
##  [75] highr_0.5.1               genefilter_1.50.0        
##  [77] S4Vectors_0.6.6           tseries_0.10-34          
##  [79] foreach_1.4.3             nortest_1.0-4            
##  [81] permute_0.8-4             boot_1.3-17              
##  [83] BiocParallel_1.2.22       bibtex_0.4.0             
##  [85] chron_2.3-47              GenomeInfoDb_1.4.3       
##  [87] matrixStats_0.15.0        moments_0.14             
##  [89] bitops_1.0-6              dmt_0.8.20               
##  [91] rgl_0.95.1367             evaluate_0.8             
##  [93] lattice_0.20-33           labeling_0.3             
##  [95] plyr_1.8.3                magrittr_1.5             
##  [97] DESeq2_1.8.2              R6_2.1.1                 
##  [99] IRanges_2.2.9             earlywarnings_1.1.22     
## [101] Hmisc_3.17-0              foreign_0.8-66           
## [103] survival_2.38-3           RCurl_1.95-4.7           
## [105] nnet_7.3-11               futile.options_1.0.0     
## [107] KernSmooth_2.23-15        ellipse_0.3-8            
## [109] locfit_1.5-9.1            grid_3.2.2               
## [111] data.table_1.9.6          Rgraphviz_2.12.0         
## [113] vegan_2.3-1               digest_0.6.8             
## [115] diptest_0.75-7            xtable_1.7-4             
## [117] stats4_3.2.2              munsell_0.4.2            
## [119] netresponse_1.18.0        quadprog_1.5-5
```

