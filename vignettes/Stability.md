### Stability analysis within group of samples

Calculate stability as the average correlation between samples and their mean for a given phylotypes vs. samples matrix:


```r
library(microbiome)

data.directory <- system.file("extdata", package = "microbiome")

genus.matrix.log10.simulated <- read.profiling(level = "oligo",		
			     	  	      method = "frpa", 
			                    data.dir = data.directory, 
			     	      	       log10 = TRUE)  

# Stability: 
stability <- estimate.stability(genus.matrix.log10.simulated)$stability
```

### Stability analysis in time

Calculate correlations between time points 1 and 2 for each subject,
then calculate stability as the average correlation of these
subject-specific correlations. Done as above, but instead of one
matrix, give two matrices as input. These are phylotypes vs. samples
matrices and correspond to time points 1 and 2, respectively. The rows
(phylotypes) and columns (subjects) should be in the same order in
both matrices.



### Version information


```r
sessionInfo()
```

```
## R version 3.1.2 (2014-10-31)
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
## [1] grid      parallel  stats     graphics  grDevices utils     datasets 
## [8] methods   base     
## 
## other attached packages:
##  [1] microbiome_0.99.34   devtools_1.7.0       googleVis_0.5.6     
##  [4] limma_3.20.9         netresponse_1.17.13  mclust_4.4          
##  [7] minet_3.20.2         infotheo_1.2.0       Rgraphviz_2.8.1     
## [10] graph_1.42.0         ggplot2_1.0.0        sorvi_0.7.13        
## [13] dplyr_0.3.0.2        rdryad_0.1.1         knitr_1.8           
## [16] gdata_2.13.3         AnnotationDbi_1.26.1 GenomeInfoDb_1.0.2  
## [19] Biobase_2.24.0       BiocGenerics_0.10.0  RSQLite_1.0.0       
## [22] DBI_0.3.1            reshape_0.8.5        vegan_2.2-1         
## [25] lattice_0.20-29      permute_0.8-3        e1071_1.6-4         
## [28] rmarkdown_0.3.10    
## 
## loaded via a namespace (and not attached):
##  [1] acepack_1.3-3.3       ape_3.1-4             assertthat_0.1       
##  [4] class_7.3-11          cluster_1.15.3        codetools_0.2-9      
##  [7] colorspace_1.2-4      df2json_0.0.2         digest_0.6.4         
## [10] dmt_0.8.20            doParallel_1.0.8      dynamicTreeCut_1.62  
## [13] evaluate_0.5.5        fastcluster_1.1.15    foreach_1.4.2        
## [16] foreign_0.8-61        formatR_1.0           Formula_1.1-2        
## [19] GO.db_2.14.0          gtable_0.1.2          gtools_3.4.1         
## [22] Hmisc_3.14-5          htmltools_0.2.6       httr_0.5             
## [25] igraph_0.7.1          impute_1.38.1         IRanges_1.22.10      
## [28] iterators_1.0.7       labeling_0.3          latticeExtra_0.6-26  
## [31] lazyeval_0.1.9        magrittr_1.0.1        MASS_7.3-37          
## [34] Matrix_1.1-4          matrixStats_0.10.3    mgcv_1.8-3           
## [37] mixOmics_5.0-3        munsell_0.4.2         mvtnorm_1.0-0        
## [40] nlme_3.1-118          nnet_7.3-8            OAIHarvester_0.1-7   
## [43] pheatmap_0.7.7        plyr_1.8.1            preprocessCore_1.26.1
## [46] proto_0.3-10          qvalue_1.38.0         RColorBrewer_1.0-5   
## [49] Rcpp_0.11.3           RCurl_1.95-4.3        reshape2_1.4.1       
## [52] RGCCA_2.0             rgl_0.95.1158         rjson_0.2.15         
## [55] RJSONIO_1.3-0         R.methodsS3_1.6.1     rpart_4.1-8          
## [58] scales_0.2.4          splines_3.1.2         stats4_3.1.2         
## [61] stringr_0.6.2         survival_2.37-7       tcltk_3.1.2          
## [64] tools_3.1.2           WGCNA_1.43            XML_3.98-1.1         
## [67] yaml_2.1.13
```

