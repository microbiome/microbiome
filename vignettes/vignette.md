---
title: "microbiome vignette"
author: "Leo Lahti and Jarkko Salojarvi"
date: "2014-09-15"
output: rmarkdown::html_vignette
vignette: > 
%\VignetteEngine{knitr::rmarkdown}
%\VignetteIndexEntry{An R Markdown Vignette made with knitr}
%\usepackage[utf8]{inputenc}
---




microbiome R package
===========

The microbiome package contains general-purpose tools for
microarray-based analysis of microbiome profiling data sets. 


### Installing and loading the release version


```r
source("http://www.bioconductor.org/biocLite.R")
biocLite("microbiome")
library(microbiome)
```

### Installing and loading the experimental development version


```r
install.packages("devtools")
```

```
## Installing package into '/home/lei/R/x86_64-pc-linux-gnu-library/3.1'
## (as 'lib' is unspecified)
```

```
## --- Please select a CRAN mirror for use in this session ---
## 
## The downloaded source packages are in
## 	'/tmp/RtmpQ4lXwS/downloaded_packages'
```

```r
library(devtools)
```

```
## 
## Attaching package: 'devtools'
## 
## The following objects are masked from 'package:utils':
## 
##     ?, help
## 
## The following object is masked from 'package:base':
## 
##     system.file
```

```r
install_github("microbiome", "microbiome")
```

```
## Installing github repo microbiome/master from microbiome
## Downloading master.zip from https://github.com/microbiome/microbiome/archive/master.zip
## Installing package from /tmp/RtmpQ4lXwS/master.zip
## arguments 'minimized' and 'invisible' are for Windows only
## Installing microbiome
## '/usr/lib/R/bin/R' --vanilla CMD INSTALL  \
##   '/tmp/RtmpQ4lXwS/devtools1520588c6168/microbiome-master'  \
##   --library='/home/lei/R/x86_64-pc-linux-gnu-library/3.1'  \
##   --install-tests
```

```r
library(microbiome)
```

```
## Loading required package: e1071
## Loading required package: vegan
## Loading required package: permute
## 
## Attaching package: 'permute'
## 
## The following object is masked from 'package:devtools':
## 
##     check
## 
## Loading required package: lattice
## This is vegan 2.0-10
## Loading required package: reshape
## 
## microbiome R package (microbiome.github.com)
## Copyright (C) 2011-2014 Leo Lahti and Jarkko Salojarvi <microbiome-admin@googlegroups.com>
## 
## 
## Attaching package: 'microbiome'
## 
## The following object is masked from 'package:vegan':
## 
##     diversity
## 
## The following object is masked from 'package:lattice':
## 
##     densityplot
## 
## The following object is masked from 'package:e1071':
## 
##     impute
```

### Examples

Further installation and usage instructions can be found at the
project [wiki](https://github.com/microbiome/microbiome/wiki/). 


### Example data set

The microbiome package contains an example data set from Lahti et al. [PeerJ 1:e32, 2013](https://peerj.com/articles/32/) concerning associations between human intestinal microbiota and blood serum lipids. Load the data in R as follows:



```r
data(peerj32)
names(peerj32)
```

```
## [1] "lipids"   "microbes" "meta"
```


### Licensing and Citations

This work can be freely used, modified and distributed under the 
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD\_licenses).

Kindly cite the work, if appropriate, as 'Leo Lahti and Jarkko Salojarvi (2014). URL: http://microbiome.github.com'. 


### References

The package utilizes tools from a number of other CRAN and
Bioconductor extensions, including ade4, df2json, rjson, fastcluster,
ggplot2, MASS, methods, minet, mixOmics, plyr, qvalue, RCurl,
reshape2, RPA, vegan, and WGCNA. We thank all authors for these
contributions:

 * N. Caballero (2013). [df2json: Convert a dataframe to JSON](http://CRAN.R-project.org/package=df2json)

 * A. Couture-Beil (2013). [rjson: JSON for R](http://CRAN.R-project.org/package=rjson)

 * A. Dabney, John D. Storey and with assistance from Gregory R. Warnes. qvalue: Q-value estimation for false discovery rate control. 

 * S. Dray and A. B. Dufour, (2007): The ade4 package: implementing the duality diagram for ecologists. Journal of Statistical Software. 22(4): 1-20.

 * S. Dejean et al. (2013). [mixOmics: Omics Data Integration Project](http://CRAN.R-project.org/package=mixOmics)

 * L. Lahti et al. A fully scalable online-preprocessing algorithm for short oligonucleotide microarray atlases. [NAR 41(10):e110, 2013](http://nar.oxfordjournals.org/content/41/10/e110)

 * L. Lahti et al. Analysis of Probe Reliability in Differential Gene Expression Studies with Short Oligonucleotide Arrays. [TCBB/IEEE 8(1):217-225, 2011](http://www.computer.org/portal/web/csdl/doi/10.1109/TCBB.2009.38)

 * L. Lahti et al. Associations between the human intestinal microbiota, Lactobacillus rhamnosus GG and serum lipids indicated by integrated analysis of high-throughput profiling data. [PeerJ 1:e32, 2013](http://dx.doi.org/10.7717/peerj.32).

 * D. T. Lang (2013). [RCurl: General network (HTTP/FTP/...) client interface for R](http://CRAN.R-project.org/package=RCurl)

 * P. Langfelder and S. Horvath, WGCNA: an R package for weighted correlation network analysis. BMC Bioinformatics 2008, 9:559 

 * P. Langfelder, S. Horvath (2012). Fast R Functions for Robust Correlations and Hierarchical Clustering. [Journal of Statistical Software, 46(11), 1-17](http://www.jstatsoft.org/v46/i11/)

 * P. E. Meyer, Frederic Lafitte and Gianluca Bontempi (2008). MINET: An open source R/Bioconductor Package for Mutual Information based Network Inference. [BMC Bioinformatics](http://www.biomedcentral.com/1471-2105/9/461)

 * D. Mullner (2013). fastcluster: Fast Hierarchical, Agglomerative Clustering Routines for R and Python. [Journal of Statistical Software, 53(9), 1-18](http://www.jstatsoft.org/v53/i09/)

 * Jari Oksanen et al. (2013). [vegan: Community Ecology Package](http://CRAN.R-project.org/package=vegan)

 * R Core Team (2013). R: A language and environment for statistical computing. R Foundation for Statistical Computing, Vienna, Austria. [ISBN 3-900051-07-0](http://www.R-project.org/)

 * W. N. Venables and B. D. Ripley (2002) Modern Applied Statistics with S. Fourth Edition. Springer, New York. ISBN 0-387-95457-0

 * H. Wickham (2007). Reshaping Data with the reshape Package. [Journal of Statistical Software, 21(12), 1-20](http://www.jstatsoft.org/v21/i12/)

 * H. Wickham. ggplot2: elegant graphics for data analysis. Springer New York, 2009.

 * H. Wickham (2011). The Split-Apply-Combine Strategy for Data Analysis. [Journal of Statistical Software, 40(1), 1-29](http://www.jstatsoft.org/v40/i01/)

### Session info


This vignette was created with


```r
sessionInfo()
```

```
## R version 3.1.1 (2014-07-10)
## Platform: x86_64-pc-linux-gnu (64-bit)
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
## [1] stats     graphics  grDevices utils     datasets  methods   base     
## 
## other attached packages:
## [1] microbiome_0.99.0 reshape_0.8.5     vegan_2.0-10      lattice_0.20-29  
## [5] permute_0.8-3     e1071_1.6-4       devtools_1.5      knitr_1.6        
## 
## loaded via a namespace (and not attached):
##  [1] acepack_1.3-3.3     class_7.3-11        cluster_1.15.3     
##  [4] codetools_0.2-9     colorspace_1.2-4    df2json_0.0.2      
##  [7] digest_0.6.4        doParallel_1.0.8    dynamicTreeCut_1.62
## [10] evaluate_0.5.5      flashClust_1.01-2   foreach_1.4.2      
## [13] foreign_0.8-61      formatR_1.0         Formula_1.1-2      
## [16] ggplot2_1.0.0       grid_3.1.1          gtable_0.1.2       
## [19] Hmisc_3.14-5        httr_0.5            igraph_0.7.1       
## [22] impute_1.38.1       iterators_1.0.7     latticeExtra_0.6-26
## [25] MASS_7.3-34         matrixStats_0.10.0  memoise_0.2.1      
## [28] minet_3.20.1        mixOmics_5.0-3      munsell_0.4.2      
## [31] nnet_7.3-8          parallel_3.1.1      pheatmap_0.7.7     
## [34] plyr_1.8.1          proto_0.3-10        RColorBrewer_1.0-5 
## [37] Rcpp_0.11.2         RCurl_1.95-4.3      reshape2_1.4       
## [40] RGCCA_2.0           rgl_0.94.1131       rjson_0.2.14       
## [43] R.methodsS3_1.6.1   rpart_4.1-8         scales_0.2.4       
## [46] splines_3.1.1       stringr_0.6.2       survival_2.37-7    
## [49] tcltk_3.1.1         tools_3.1.1         WGCNA_1.41-1       
## [52] whisker_0.4
```




