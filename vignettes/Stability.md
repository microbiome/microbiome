
### Stability analysis 

Calculate stability as the average correlation between samples and their
mean for a given phylotypes vs. samples matrix.

Download example data (Lahti et al. Nat. Comm. 5:4344, 2014):


```r
library(microbiome)
pseq <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

Quantify temporal stability across the abundance range for each taxa


```r
stability <- intermediate_stability(pseq, method = "correlation")
intermediate.stability <- sapply(stability, function (x) {x$stability})
```


## Bimodality analysis

Multiple scores are available to calculate coefficients of bimodality
for each taxonomic group. Some of the scores are more generally
quantifying multimodality, or deviation from unimodality but do not
yield estimates of the number of modes. 


```r
# Pick samples from zero time point (cross-sectional analysis)
pseq0 <- subset_samples(pseq, time == 0 & DNA_extraction_method == "r")

# Bimodality score with potential analysis + bootstrap 
bimodality <- multimodality(pseq0, method = "potential_bootstrap")

# Coefficient of bimodality. Also see the function coefficient_of_bimodality
bimodality.coef <- multimodality(pseq0, method = "coefficient_of_bimodality")
```

DIP test for multimodality:


```r
# DIP multimodality test uses the separate diptest package.
library(diptest)

# Pick OTU log10 data
otu <- otu_table(pseq0)@.Data
dip <- apply(otu, 1, function (x) dip.test(x, simulate.p.value = TRUE, B = 1000))
dip2 <- data.frame(t(sapply(dip, function (x) {c(x$statistic, x$p.value)})))
colnames(dip2) <- c("score", "p.value")
dip2$tax <- names(dip)

# Dip measures unimodality. Values range between 0 to 1. 
# Values less than 0.05 indicate significant deviation from unimodality. 
# To score multimodality, use the inverse:
multimodality.dip <- 1 - dip2$score
```


Compare the alternative bimodality scores


```r
pairs(cbind(pb = bimodality, bc = bimodality.coef, dip = multimodality.dip))
```

```
## Error in pairs.default(cbind(pb = bimodality, bc = bimodality.coef, dip = multimodality.dip)): non-numeric argument to 'pairs'
```


Visualize population densities for selected taxa


```r
# Pick the most and least bimodal taxa as examples
unimodal <- names(which.min(bimodality))
bimodal <- names(which.max(bimodality))

# Visualize population frequencies
p1 <- plot_density(pseq, otu.name = unimodal, log10 = TRUE) 
p2 <- plot_density(pseq, otu.name = bimodal, log10 = TRUE) 
library(gridExtra)
grid.arrange(p1, p2, nrow = 1)
```

![plot of chunk stability2](figure/stability2-1.png) 


## Bimodality versus intermediate stability

The analysis suggests that bimodal taxa tend to have instable intermediate abundances (following [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) but with a different subset, parameters and model):


```r
# For clarity, visualize only the prevalent taxa seen with HITChip signal >250
# (note HITChip signal reflects read count but is conceptually different) 
# in at least 10% of the samples. 
pseqf <- filter_taxa(pseq0, function(x) sum(x > 300) > (0.2*length(x)), TRUE)

# Taxa that have bistability and bimodality estimates
s <- taxa_names(pseqf)

plot(intermediate.stability[s], bimodality[s], xlab = "Intermediate stability", ylab = "Bimodality", type = "n")
```

```
## Error in plot.window(...): need finite 'ylim' values
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png) 

```r
text(intermediate.stability[s], bimodality[s], label = s, cex = 0.7)
```


## Variability within group of samples (inter-individual stability)

Assess 'inter-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between samples and their mean for a given samples vs phylotypes matrix. For the illustration, calculate inter-individual stability (variability) separately for Placebo and LGG groups.

Load example data


```r
library(microbiome)
x <- download_microbiome("dietswap")
# Add time field (two time points needed within each group for the 
# intraindividual method)
sample_data(x)$time <- sample_data(x)$timepoint.within.group
```


### Inter-individual variability

Variability across subjects within a group


```r
res <- estimate_variability(x, "interindividual")
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Inter-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

![plot of chunk variability-example2d](figure/variability-example2d-1.png) 


### Intra-individual variability

Variability within subjects over time (intra-individual stability). Assess 'intra-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between two time points within subjects, calculated separately within each group. For illustration, check intra-individual stability (variability) separately for Placebo and LGG groups.


```r
res <- estimate_variability(x, "intraindividual")
```


Visualize


```r
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Intra-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

```
## Warning in loop_apply(n, do.ply): Removed 6 rows containing non-finite
## values (stat_boxplot).
```

![plot of chunk variability-intra](figure/variability-intra-1.png) 


### TODO

Add examples on tipping elements.

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
## [1] grid      stats     graphics  grDevices utils     datasets  methods  
## [8] base     
## 
## other attached packages:
##  [1] gridExtra_0.9.1    ggplot2_1.0.1      diptest_0.75-6    
##  [4] dplyr_0.4.1        limma_3.24.3       knitr_1.10.5      
##  [7] ade4_1.7-2         phyloseq_1.13.2    microbiome_0.99.48
## [10] vegan_2.2-1        lattice_0.20-31    permute_0.8-3     
## [13] scimapClient_0.2.1
## 
## loaded via a namespace (and not attached):
##  [1] nlme_3.1-120              bitops_1.0-6             
##  [3] RColorBrewer_1.1-2        GenomeInfoDb_1.4.0       
##  [5] tools_3.2.0               rpart_4.1-9              
##  [7] KernSmooth_2.23-14        lazyeval_0.1.10          
##  [9] nortest_1.0-3             Hmisc_3.16-0             
## [11] DBI_0.3.1                 BiocGenerics_0.14.0      
## [13] mgcv_1.8-6                colorspace_1.2-6         
## [15] nnet_7.3-9                DESeq2_1.8.1             
## [17] moments_0.14              chron_2.3-45             
## [19] rdryad_0.1.1              Biobase_2.28.0           
## [21] formatR_1.2               RGCCA_2.0                
## [23] labeling_0.3              tseries_0.10-34          
## [25] scales_0.2.4              lmtest_0.9-33            
## [27] genefilter_1.50.0         quadprog_1.5-5           
## [29] OAIHarvester_0.1-7        tgp_2.4-11               
## [31] stringr_1.0.0             digest_0.6.8             
## [33] foreign_0.8-63            earlywarnings_1.1.19     
## [35] XVector_0.8.0             maps_2.3-9               
## [37] RSQLite_1.0.0             zoo_1.7-12               
## [39] gtools_3.4.2              BiocParallel_1.2.1       
## [41] acepack_1.3-3.3           RCurl_1.95-4.6           
## [43] magrittr_1.5              Formula_1.2-1            
## [45] futile.logger_1.4.1       Matrix_1.2-0             
## [47] Rcpp_0.11.6               munsell_0.4.2            
## [49] S4Vectors_0.6.0           maptree_1.4-7            
## [51] ape_3.2                   proto_0.3-10             
## [53] stringi_0.4-1             MASS_7.3-40              
## [55] RJSONIO_1.3-0             zlibbioc_1.14.0          
## [57] plyr_1.8.2                gdata_2.16.1             
## [59] parallel_3.2.0            Biostrings_2.36.1        
## [61] splines_3.2.0             multtest_2.24.0          
## [63] annotate_1.46.0           locfit_1.5-9.1           
## [65] igraph_0.7.1              GenomicRanges_1.20.3     
## [67] boot_1.3-16               geneplotter_1.46.0       
## [69] reshape2_1.4.1            codetools_0.2-11         
## [71] mixOmics_5.0-4            stats4_3.2.0             
## [73] futile.options_1.0.0      XML_3.98-1.1             
## [75] evaluate_0.7              RcppArmadillo_0.5.100.1.0
## [77] latticeExtra_0.6-26       biom_0.3.12              
## [79] lambda.r_1.1.7            data.table_1.9.4         
## [81] spam_1.0-1                foreach_1.4.2            
## [83] gtable_0.1.2              assertthat_0.1           
## [85] xtable_1.7-4              e1071_1.6-4              
## [87] sorvi_0.7.23              class_7.3-12             
## [89] Kendall_2.2               survival_2.38-1          
## [91] pheatmap_1.0.2            iterators_1.0.7          
## [93] som_0.3-5                 AnnotationDbi_1.30.1     
## [95] IRanges_2.2.1             fields_8.2-1             
## [97] cluster_2.0.1             rgl_0.95.1247
```

