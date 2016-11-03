---
title: "Stability"
author: "Leo Lahti"
date: "2016-11-03"
bibliography: 
- bibliography.bib
- references.bib
output: 
  rmarkdown::html_vignette
---
<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - stability}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->


# Microbiome stability analysis

Get example data - [HITChip Atlas of 130 genus-like taxa across 1006 healthy western adults](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). A subset of 76 subjects have also short time series available for temporal stability analysis:


```r
library(microbiome)
data(atlas1006)
pseq <- atlas1006

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
al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)
used straightforward correlation analysis to quantify how the distance
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

Some other standard multimodality tests include the DIP test from the
[diptest](https://cran.r-project.org/web/packages/diptest/index.html)
package.



Visualize population densities 


```r
# Pick the most and least bimodal taxa as examples
bimodality <- bimodality.pb
unimodal <- names(which.min(bimodality))
bimodal  <- names(which.max(bimodality))

# Visualize population frequencies
library(ggplot2)
theme_set(theme_bw(20))
p1 <- plot_density(pseq, variable = unimodal, log10 = TRUE) 
p2 <- plot_density(pseq, variable = bimodal,  log10 = TRUE) 
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

# Repel overlapping labels
# Install ggrepel package if needed
# install.packages("devtools")
# devtools::install_github("slowkow/ggrepel")
# See https://github.com/slowkow/ggrepel/blob/master/vignettes/ggrepel.md
library(ggrepel)
p <- p + geom_text_repel(size = 3)

print(p)
```

![plot of chunk bimodalitybistability](figure/bimodalitybistability-1.png)

## Detecting a tipping point

Identify potential minima in cross-section population data as
tipping point candidates. 


```r
# Pick example data
library(phyloseq)
library(microbiome)
data("atlas1006")
pseq <- atlas1006
pseq <- subset_samples(pseq, DNA_extraction_method == "r")
pseq <- transform_phyloseq(pseq, "relative.abundance")

# Dialister log10 relative abundance
x <- log10(get_sample(pseq, "Dialister"))

# Potential analysis to identify potential minima
library(earlywarnings)
res <- livpotential_ews(x)

# Identify the potential minimum location as a tipping point candidate 
tipping.point <- 10^res$min.points

print(tipping.point)
```

```
## [1] 0.2491531
```

## Variation lineplot and Bimodality hotplot

Pick subset of the [HITChip Atlas data set](http://doi.org/10.5061/dryad.pk75d) and plot the subject abundance variation lineplot (**Variation lineplot**) and **Bimodality hotplot** for a given taxon as in [Lahti et al. 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). The bi-stable Dialister has bimodal population distribution and reduced temporal stability within subjects at intermediate abundances.


```r
# Variation line plot:
# Indicates the abundance variation range
# for subjects with multiple time points
pv <- plot_variation(pseq, "Dialister", tipping.point = tipping.point, xlim = c(0.01, 100))
print(pv)

# Bimodality hotplot:
# Only consider a unique sample from each subject
# baseline time point for density plot
pseq.baseline <- subset_samples(pseq, time == 0)
ph <- plot_bimodal(pseq.baseline, "Dialister", tipping.point = tipping.point)
print(ph)
```

<img src="figure/stability-variationplot-1.png" title="plot of chunk stability-variationplot" alt="plot of chunk stability-variationplot" width="430px" /><img src="figure/stability-variationplot-2.png" title="plot of chunk stability-variationplot" alt="plot of chunk stability-variationplot" width="430px" />



### Potential analysis

Potential analysis, used for instance in [Hirota et al. Science, 334, 232-235.](http://www.sciencemag.org/content/334/6053/232.long)


```r
# Create simulated example data
X <- c(rnorm(1000, mean = 0), rnorm(1000, mean = -2), 
 	           rnorm(1000, mean = 2))
param <- seq(0,5,length=3000) 

# Run potential analysis
res <- movpotential(X, param)
```

```
## Error in eval(expr, envir, enclos): could not find function "movpotential"
```

```r
# Visualize
p <- plot_potential(res$res, title = '', 
	       	   xlab.text = '', ylab.text = '', 
		   cutoff = 0.5, plot.contours = TRUE, binwidth = 0.2)
```

```
## Error in plot_potential(res$res, title = "", xlab.text = "", ylab.text = "", : unused arguments (title = "", xlab.text = "", ylab.text = "")
```

```r
print(p)
```

![plot of chunk movpotential](figure/movpotential-1.png)



### Version information


```r
sessionInfo()
```

```
## R version 3.3.2 (2016-10-31)
## Platform: x86_64-pc-linux-gnu (64-bit)
## Running under: Ubuntu 16.04.1 LTS
## 
## locale:
##  [1] LC_CTYPE=en_US.UTF-8       LC_NUMERIC=C              
##  [3] LC_TIME=de_BE.UTF-8        LC_COLLATE=en_US.UTF-8    
##  [5] LC_MONETARY=de_BE.UTF-8    LC_MESSAGES=en_US.UTF-8   
##  [7] LC_PAPER=de_BE.UTF-8       LC_NAME=C                 
##  [9] LC_ADDRESS=C               LC_TELEPHONE=C            
## [11] LC_MEASUREMENT=de_BE.UTF-8 LC_IDENTIFICATION=C       
## 
## attached base packages:
##  [1] tcltk     splines   stats4    grid      parallel  stats     graphics 
##  [8] grDevices utils     datasets  methods   base     
## 
## other attached packages:
##  [1] earlywarnings_1.1.22  tseries_0.10-35       tgp_2.4-14           
##  [4] moments_0.14          ggrepel_0.6.3         gridExtra_2.2.1      
##  [7] HITChipDB_0.6.30      RPA_1.28.0            affy_1.50.0          
## [10] Biobase_2.32.0        RMySQL_0.10.9         preprocessCore_1.34.0
## [13] DBI_0.5-1             scales_0.4.0          SpiecEasi_0.1        
## [16] VGAM_1.0-2            huge_1.2.7            igraph_1.0.1         
## [19] Matrix_1.2-7.1        FD_1.0-12             vegan_2.4-1          
## [22] lattice_0.20-34       permute_0.9-4         geometry_0.3-6       
## [25] magic_1.5-6           abind_1.4-5           ape_3.5              
## [28] ade4_1.7-4            RColorBrewer_1.1-2    knitcitations_1.0.7.1
## [31] knitr_1.14            intergraph_2.0-2      sna_2.4              
## [34] statnet.common_3.3.0  network_1.13.0        ggnet_0.1.0          
## [37] GGally_1.2.0          devtools_1.12.0       limma_3.28.21        
## [40] sorvi_0.7.46          tibble_1.2            ggplot2_2.1.0        
## [43] tidyr_0.6.0           dplyr_0.5.0           MASS_7.3-45          
## [46] netresponse_1.32.2    reshape2_1.4.1        mclust_5.2           
## [49] minet_3.30.0          Rgraphviz_2.16.0      graph_1.50.0         
## [52] BiocGenerics_0.18.0   microbiome_0.99.87    phyloseq_1.16.2      
## 
## loaded via a namespace (and not attached):
##  [1] colorspace_1.2-7      dynamicTreeCut_1.63-1 som_0.3-5.1          
##  [4] qvalue_2.4.2          XVector_0.12.1        affyio_1.42.0        
##  [7] AnnotationDbi_1.34.4  mvtnorm_1.0-5         lubridate_1.6.0      
## [10] RefManageR_0.11.0     codetools_0.2-15      doParallel_1.0.10    
## [13] impute_1.46.0         spam_1.4-0            Formula_1.2-1        
## [16] jsonlite_1.1          Cairo_1.5-9           WGCNA_1.51           
## [19] cluster_2.0.5         GO.db_3.3.0           Kendall_2.2          
## [22] httr_1.2.1            assertthat_0.1        lazyeval_0.2.0       
## [25] formatR_1.4           acepack_1.3-3.3       tools_3.3.2          
## [28] gtable_0.2.0          maps_3.1.1            Rcpp_0.12.7          
## [31] Biostrings_2.40.2     RJSONIO_1.3-0         multtest_2.28.0      
## [34] nlme_3.1-128          iterators_1.0.8       lmtest_0.9-34        
## [37] fastcluster_1.1.21    stringr_1.1.0         XML_3.98-1.4         
## [40] zoo_1.7-13            zlibbioc_1.18.0       BiocInstaller_1.22.3 
## [43] biomformat_1.0.2      rhdf5_2.16.0          fields_8.4-1         
## [46] memoise_1.0.0         rpart_4.1-10          reshape_0.8.5        
## [49] latticeExtra_0.6-28   stringi_1.1.2         maptree_1.4-7        
## [52] RSQLite_1.0.0         highr_0.6             S4Vectors_0.10.3     
## [55] nortest_1.0-4         foreach_1.4.3         boot_1.3-18          
## [58] bibtex_0.4.0          chron_2.3-47          matrixStats_0.51.0   
## [61] bitops_1.0-6          dmt_0.8.20            evaluate_0.10        
## [64] labeling_0.3          plyr_1.8.4            magrittr_1.5         
## [67] R6_2.2.0              IRanges_2.6.1         Hmisc_3.17-4         
## [70] foreign_0.8-67        withr_1.0.2           mgcv_1.8-15          
## [73] survival_2.39-5       RCurl_1.95-4.8        nnet_7.3-12          
## [76] KernSmooth_2.23-15    data.table_1.9.6      digest_0.6.10        
## [79] munsell_0.4.3         quadprog_1.5-5
```

