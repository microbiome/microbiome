
### Stability analysis 

Calculate stability as the average correlation between samples and their
mean for a given phylotypes vs. samples matrix:


```r
# Download example data (Lahti et al. Nat. Comm. 5:4344, 2014)
library(microbiome)
data.atlas1006 <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```r
# Quantify temporal stability across the abundance range
stability <- intermediate_stability(data.atlas1006, method = "correlation")
intermediate.stability <- sapply(stability, function (x) {x$stability})
hist(intermediate.stability, main = "Intermediate stability")
```

![plot of chunk bistability](figure/bistability-1.png) 


## Bimodality analysis

Multiple scores are available to calculate coefficients of bimodality
for each taxonomic group. Some of the scores are more generally
quantifying multimodality, or deviation from unimodality but do not
yield estimates of the number of modes. 


```r
# Pick samples from zero time point (cross-sectional analysis)
pseq0 <- subset_samples(data.atlas1006, time == 0 & DNA_extraction_method == "r")

# Pick OTU log10 data
otu <- otu_table(pseq0)@.Data
otu.log10 <- log10(otu)

# Calculate coefficient of bimodality for each taxa (calculated for log10 abundances)
# Potential analysis + Bootstrap bimodality score
bimodality <- multimodality(pseq0, method = "potential.bootstrap")

# Coefficient of bimodality. Also see the function coefficient_of_bimodality
bimodality.coef <- multimodality(pseq0, method = "coefficient_of_bimodality")

# Alo calculate DIP multimodality test. This uses the separate diptest package.
# This seems to coincide much better with the other scores when absolute OTU count is used
dip <- apply(otu, 1, function (x) dip.test(x, simulate.p.value = TRUE, B = 5000))
```

```
## Error in FUN(newX[, i], ...): could not find function "dip.test"
```

```r
dip2 <- data.frame(t(sapply(dip, function (x) {c(x$statistic, x$p.value)})))
```

```
## Error in t(sapply(dip, function(x) {: error in evaluating the argument 'x' in selecting a method for function 't': Error in sapply(dip, function(x) { : 
##   error in evaluating the argument 'X' in selecting a method for function 'sapply': Error: object 'dip' not found
```

```r
dip2$tax <- names(dip)
```

```
## Error in eval(expr, envir, enclos): object 'dip' not found
```

```r
colnames(dip2) <- c("score", "p.value")
```

```
## Error in colnames(dip2) <- c("score", "p.value"): object 'dip2' not found
```

```r
# Dip measures unimodality. Values range between 0 to 1. Values less than 0.05 indicate significant 
# deviation from unimodality. To score multimodality, use the inverse:
multimodality.dip <- 1 - dip2$score
```

```
## Error in eval(expr, envir, enclos): object 'dip2' not found
```

```r
#multimodality.dip <- dip2$p.value

# Compare the alternative bimodality scores
pairs(cbind(pb = bimodality, bc = bimodality.coef, dip = multimodality.dip))
```

```
## Error in eval(expr, envir, enclos): object 'multimodality.dip' not found
```

Visualize population densities for selected taxa


```r
# Pick the most and least bimodal taxa as examples
unimodal <- names(which.min(bimodality))
bimodal <- names(which.max(bimodality))

# Visualize population frequencies
par(mfrow = c(1,2))
plot(density(otu.log10[unimodal, ]), main = unimodal)
plot(density(otu.log10[bimodal, ]), main = bimodal)
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
text(intermediate.stability[s], bimodality[s], label = s, cex = 0.7)
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png) 


### Variability within group of samples (inter-individual stability)

Assess 'inter-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between samples and their mean for a given samples vs phylotypes matrix. For the illustration, calculate inter-individual stability (variability) separately for Placebo and LGG groups.


```r
# Example data
library(microbiome)
x <- download_microbiome("dietswap")

# Add time field (two time points needed within each group for the 
# intraindividual method)
sample_data(x)$time <- sample_data(x)$timepoint.within.group

# Estimate inter-individual variability
res <- estimate_variability(x, "interindividual")

# Estimate intra-individual variability
res <- estimate_variability(x, "intraindividual")

# Visualize
library(ggplot2)
theme_set(theme_bw(20))
p <- ggplot(res$data, aes(x = group, y = correlation))
p <- p + geom_boxplot()
p <- p + ggtitle(paste("Inter-individual variability (p=", round(res$p.value, 6), ")"))
p <- p + ylab("Correlation")
print(p)
```

```
## Warning in loop_apply(n, do.ply): Removed 6 rows containing non-finite
## values (stat_boxplot).
```

![plot of chunk variability-example2](figure/variability-example2-1.png) 


### Variability within subjects over time (intra-individual stability)

Assess 'intra-individual stability' as in [Salonen et al. ISME J 2014](http://www.nature.com/ismej/journal/v8/n11/full/ismej201463a.html). This is defined as the average correlation between two time points within subjects, calculated separately within each group. For illustration, check intra-individual stability (variability) separately for Placebo and LGG groups.


```r
res <- estimate_variability(x, m, "intraindividual")
```

```
## Warning in if (type == "interindividual") {: the condition has length > 1
## and only the first element will be used
```

```
## Warning in if (type == "intraindividual") {: the condition has length > 1
## and only the first element will be used
```

```
## Error in estimate_variability(x, m, "intraindividual"): object 'variability' not found
```

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
##  [1] tcltk     parallel  grid      stats     graphics  grDevices utils    
##  [8] datasets  methods   base     
## 
## other attached packages:
##  [1] googleVis_0.5.8       HITChipDB_0.5.16      RPA_1.24.0           
##  [4] affy_1.46.0           Biobase_2.28.0        BiocGenerics_0.14.0  
##  [7] RMySQL_0.10.3         preprocessCore_1.30.0 ade4_1.7-2           
## [10] limma_3.24.3          scales_0.2.4          RSQLite_1.0.0        
## [13] DBI_0.3.1             netresponse_1.18.0    reshape_0.8.5        
## [16] mclust_5.0.1          minet_3.26.0          Rgraphviz_2.12.0     
## [19] graph_1.46.0          ggplot2_1.0.1         sorvi_0.7.23         
## [22] reshape2_1.4.1        dplyr_0.4.1           phyloseq_1.13.2      
## [25] microbiome_0.99.48    vegan_2.2-1           lattice_0.20-31      
## [28] permute_0.8-3         rdryad_0.1.1          knitcitations_1.0.5  
## [31] knitr_1.10.5          rmarkdown_0.6.1       scimapClient_0.2.1   
## 
## loaded via a namespace (and not attached):
##   [1] spam_1.0-1                Hmisc_3.16-0             
##   [3] plyr_1.8.2                igraph_0.7.1             
##   [5] lazyeval_0.1.10           splines_3.2.0            
##   [7] BiocParallel_1.2.1        GenomeInfoDb_1.4.0       
##   [9] digest_0.6.8              foreach_1.4.2            
##  [11] BiocInstaller_1.18.2      htmltools_0.2.6          
##  [13] GO.db_3.1.2               gdata_2.16.1             
##  [15] magrittr_1.5              memoise_0.2.1            
##  [17] cluster_2.0.1             doParallel_1.0.8         
##  [19] fastcluster_1.1.16        Biostrings_2.36.1        
##  [21] annotate_1.46.0           matrixStats_0.14.0       
##  [23] Kendall_2.2               tseries_0.10-34          
##  [25] colorspace_1.2-6          RCurl_1.95-4.6           
##  [27] RcppArmadillo_0.5.100.1.0 genefilter_1.50.0        
##  [29] impute_1.42.0             survival_2.38-1          
##  [31] zoo_1.7-12                iterators_1.0.7          
##  [33] ape_3.2                   gtable_0.1.2             
##  [35] zlibbioc_1.14.0           XVector_0.8.0            
##  [37] RGCCA_2.0                 tgp_2.4-11               
##  [39] maps_2.3-9                futile.options_1.0.0     
##  [41] pheatmap_1.0.2            mvtnorm_1.0-2            
##  [43] som_0.3-5                 bibtex_0.4.0             
##  [45] Rcpp_0.11.6               xtable_1.7-4             
##  [47] foreign_0.8-63            Formula_1.2-1            
##  [49] stats4_3.2.0              httr_0.6.1               
##  [51] RColorBrewer_1.1-2        acepack_1.3-3.3          
##  [53] XML_3.98-1.1              nnet_7.3-9               
##  [55] locfit_1.5-9.1            RJSONIO_1.3-0            
##  [57] dynamicTreeCut_1.62       labeling_0.3             
##  [59] AnnotationDbi_1.30.1      munsell_0.4.2            
##  [61] tools_3.2.0               moments_0.14             
##  [63] evaluate_0.7              stringr_1.0.0            
##  [65] dmt_0.8.20                maptree_1.4-7            
##  [67] RefManageR_0.8.45         rgl_0.95.1247            
##  [69] nlme_3.1-120              formatR_1.2              
##  [71] e1071_1.6-4               affyio_1.36.0            
##  [73] geneplotter_1.46.0        stringi_0.4-1            
##  [75] futile.logger_1.4.1       fields_8.2-1             
##  [77] earlywarnings_1.1.19      Matrix_1.2-0             
##  [79] multtest_2.24.0           biom_0.3.12              
##  [81] lmtest_0.9-33             data.table_1.9.4         
##  [83] bitops_1.0-6              GenomicRanges_1.20.3     
##  [85] qvalue_2.0.0              R6_2.0.1                 
##  [87] latticeExtra_0.6-26       KernSmooth_2.23-14       
##  [89] gridExtra_0.9.1           IRanges_2.2.1            
##  [91] codetools_0.2-11          lambda.r_1.1.7           
##  [93] boot_1.3-16               MASS_7.3-40              
##  [95] gtools_3.4.2              assertthat_0.1           
##  [97] chron_2.3-45              proto_0.3-10             
##  [99] DESeq2_1.8.1              OAIHarvester_0.1-7       
## [101] nortest_1.0-3             S4Vectors_0.6.0          
## [103] mgcv_1.8-6                mixOmics_5.0-4           
## [105] quadprog_1.5-5            rpart_4.1-9              
## [107] class_7.3-12              WGCNA_1.46               
## [109] lubridate_1.3.3
```

