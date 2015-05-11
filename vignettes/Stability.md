
### Stability analysis 

Calculate stability as the average correlation between samples and their
mean for a given phylotypes vs. samples matrix:


```r
# Example data
library(microbiome)
data.atlas1006 <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```r
x <- log10(data.atlas1006$microbes)
m <- data.atlas1006$meta

# Quantify temporal stability across the abundance range
stability <- estimate_stability(x, m, method = "correlation")
intermediate.stability <- sapply(stability, function (x) {x$stability})
hist(intermediate.stability, main = "Intermediate stability")
```

![plot of chunk bistability](figure/bistability-1.png) 


## Bimodality quantification

Calculate coefficient of bimodality (used in [Shade et
al.](http://mbio.asm.org/content/5/4/e01371-14)) for taxa in an
example data set, and plot the taxa with the lowest and highest score:


```r
# Pick samples from zero time point (cross-sectional analysis)
s <- subset(m, time == 0)$sample
m0 <- m[s,]
x0 <- x[s,]

# Calculate coefficient of bimodality for each taxa
bimodality <- apply(x0, 2, coefficient.of.bimodality)

# Pick the most and least bimodal taxa as examples
unimodal <- names(which.min(bimodality))
bimodal <- names(which.max(bimodality))

# Visualize
par(mfrow = c(1,2))
plot(density(x[, unimodal]), main = unimodal)
plot(density(x[, bimodal]), main = bimodal)
```

![plot of chunk bimodality](figure/bimodality-1.png) 

## Bimodality versus intermediate stability

The analysis suggests that bimodal taxa tend to have instable intermediate abundances:


```r
s <- intersect(names(intermediate.stability), names(bimodality))
plot(intermediate.stability[s], bimodality[s], xlab = "Intermediate stability", ylab = "Bimodality")
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png) 


### Variability within group of samples (inter-individual stability)

Assess 'inter-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between samples and their mean for a given samples vs phylotypes matrix. For the illustration, calculate inter-individual stability (variability) separately for Placebo and LGG groups.


```r
# Example data
library(microbiome)
data.peerj32 <- download_microbiome("peerj32")
x <- data.peerj32$microbes
m <- data.peerj32$meta

# Estimate inter-individual variability
res <- estimate_variability(x, m, "interindividual")
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Inter-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

![plot of chunk variability-example2](figure/variability-example2-1.png) 


### Variability within subjects over time (intra-individual stability)

Assess 'intra-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between two time points within subjects, calculated separately within each group. For illustration, check intra-individual stability (variability) separately for Placebo and LGG groups.


```r
res <- estimate_variability(x, m, "intraindividual")
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Intra-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

![plot of chunk variability-intra](figure/variability-intra-1.png) 



### Version information


```r
sessionInfo()
```

```
## R version 3.2.0 (2015-04-16)
## Platform: x86_64-unknown-linux-gnu (64-bit)
## Running under: Ubuntu 14.10
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
##  [1] tcltk     parallel  grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] scales_0.2.4          googleVis_0.5.8       HITChipDB_0.5.15     
##  [4] RPA_1.24.0            affy_1.46.0           Biobase_2.28.0       
##  [7] BiocGenerics_0.14.0   RMySQL_0.10.3         DBI_0.3.1            
## [10] preprocessCore_1.30.0 limma_3.24.3          gdata_2.13.3         
## [13] netresponse_1.18.0    reshape_0.8.5         mclust_5.0.1         
## [16] minet_3.26.0          Rgraphviz_2.12.0      graph_1.46.0         
## [19] ggplot2_1.0.1         sorvi_0.7.23          microbiome_0.99.46   
## [22] vegan_2.2-1           lattice_0.20-31       permute_0.8-3        
## [25] reshape2_1.4.1        phyloseq_1.12.2       e1071_1.6-4          
## [28] dplyr_0.4.1           ade4_1.7-2            rdryad_0.1.1         
## [31] knitcitations_1.0.5   knitr_1.10            rmarkdown_0.5.3.5    
## [34] scimapClient_0.2.1   
## 
## loaded via a namespace (and not attached):
##   [1] spam_1.0-1              Hmisc_3.16-0           
##   [3] plyr_1.8.2              igraph_0.7.1           
##   [5] df2json_0.0.2           lazyeval_0.1.10        
##   [7] splines_3.2.0           BiocParallel_1.2.1     
##   [9] GenomeInfoDb_1.4.0      digest_0.6.8           
##  [11] foreach_1.4.2           BiocInstaller_1.18.1   
##  [13] htmltools_0.2.6         magrittr_1.5           
##  [15] memoise_0.2.1           cluster_2.0.1          
##  [17] fastcluster_1.1.16      Biostrings_2.36.0      
##  [19] annotate_1.46.0         Kendall_2.2            
##  [21] tseries_0.10-34         colorspace_1.2-6       
##  [23] RCurl_1.95-4.6          RcppArmadillo_0.5.000.0
##  [25] genefilter_1.50.0       survival_2.38-1        
##  [27] zoo_1.7-12              iterators_1.0.7        
##  [29] ape_3.2                 gtable_0.1.2           
##  [31] zlibbioc_1.14.0         XVector_0.8.0          
##  [33] RGCCA_2.0               tgp_2.4-11             
##  [35] maps_2.3-9              futile.options_1.0.0   
##  [37] pheatmap_1.0.2          mvtnorm_1.0-2          
##  [39] som_0.3-5               bibtex_0.4.0           
##  [41] Rcpp_0.11.5             xtable_1.7-4           
##  [43] foreign_0.8-63          Formula_1.2-1          
##  [45] stats4_3.2.0            httr_0.6.1             
##  [47] RColorBrewer_1.1-2      acepack_1.3-3.3        
##  [49] XML_3.98-1.1            nnet_7.3-9             
##  [51] locfit_1.5-9.1          RJSONIO_1.3-0          
##  [53] labeling_0.3            AnnotationDbi_1.30.1   
##  [55] munsell_0.4.2           tools_3.2.0            
##  [57] moments_0.14            RSQLite_1.0.0          
##  [59] devtools_1.7.0          evaluate_0.7           
##  [61] stringr_1.0.0           dmt_0.8.20             
##  [63] maptree_1.4-7           RefManageR_0.8.45      
##  [65] rgl_0.95.1201           nlme_3.1-120           
##  [67] formatR_1.2             affyio_1.36.0          
##  [69] geneplotter_1.46.0      stringi_0.4-1          
##  [71] futile.logger_1.4.1     fields_8.2-1           
##  [73] earlywarnings_1.1.19    Matrix_1.2-0           
##  [75] multtest_2.24.0         biom_0.3.12            
##  [77] lmtest_0.9-33           data.table_1.9.4       
##  [79] bitops_1.0-6            GenomicRanges_1.20.3   
##  [81] qvalue_2.0.0            latticeExtra_0.6-26    
##  [83] KernSmooth_2.23-14      gridExtra_0.9.1        
##  [85] IRanges_2.2.1           codetools_0.2-11       
##  [87] lambda.r_1.1.7          boot_1.3-16            
##  [89] MASS_7.3-40             gtools_3.4.2           
##  [91] assertthat_0.1          chron_2.3-45           
##  [93] proto_0.3-10            DESeq2_1.8.0           
##  [95] OAIHarvester_0.1-7      rjson_0.2.15           
##  [97] nortest_1.0-3           S4Vectors_0.6.0        
##  [99] mgcv_1.8-6              mixOmics_5.0-4         
## [101] quadprog_1.5-5          rpart_4.1-9            
## [103] class_7.3-12            lubridate_1.3.3
```

