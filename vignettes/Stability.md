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
##  [1] googleVis_0.5.6      limma_3.20.9         sorvi_0.7.13        
##  [4] dplyr_0.3.0.2        gdata_2.13.3         netresponse_1.17.13 
##  [7] mclust_4.4           minet_3.20.2         infotheo_1.2.0      
## [10] Rgraphviz_2.8.1      graph_1.42.0         ggplot2_1.0.0       
## [13] microbiome_0.99.35   AnnotationDbi_1.26.1 GenomeInfoDb_1.0.2  
## [16] Biobase_2.24.0       BiocGenerics_0.10.0  RSQLite_1.0.0       
## [19] DBI_0.3.1            reshape_0.8.5        vegan_2.2-1         
## [22] lattice_0.20-29      permute_0.8-3        e1071_1.6-4         
## [25] knitr_1.8            rdryad_0.1.1         knitcitations_1.0.5 
## [28] rmarkdown_0.3.10    
## 
## loaded via a namespace (and not attached):
##  [1] acepack_1.3-3.3       ape_3.1-4             assertthat_0.1       
##  [4] bibtex_0.4.0          class_7.3-11          cluster_1.15.3       
##  [7] codetools_0.2-9       colorspace_1.2-4      df2json_0.0.2        
## [10] digest_0.6.4          dmt_0.8.20            doParallel_1.0.8     
## [13] dynamicTreeCut_1.62   evaluate_0.5.5        fastcluster_1.1.15   
## [16] foreach_1.4.2         foreign_0.8-61        formatR_1.0          
## [19] Formula_1.1-2         GO.db_2.14.0          gtable_0.1.2         
## [22] gtools_3.4.1          Hmisc_3.14-5          htmltools_0.2.6      
## [25] httr_0.5              igraph_0.7.1          impute_1.38.1        
## [28] IRanges_1.22.10       iterators_1.0.7       labeling_0.3         
## [31] latticeExtra_0.6-26   lazyeval_0.1.9        lubridate_1.3.3      
## [34] magrittr_1.0.1        MASS_7.3-37           Matrix_1.1-4         
## [37] matrixStats_0.10.3    memoise_0.2.1         mgcv_1.8-3           
## [40] mixOmics_5.0-3        munsell_0.4.2         mvtnorm_1.0-0        
## [43] nlme_3.1-118          nnet_7.3-8            OAIHarvester_0.1-7   
## [46] pheatmap_0.7.7        plyr_1.8.1            preprocessCore_1.26.1
## [49] proto_0.3-10          qvalue_1.38.0         RColorBrewer_1.0-5   
## [52] Rcpp_0.11.3           RCurl_1.95-4.3        RefManageR_0.8.45    
## [55] reshape2_1.4.1        RGCCA_2.0             rgl_0.95.1158        
## [58] rjson_0.2.15          RJSONIO_1.3-0         R.methodsS3_1.6.1    
## [61] rpart_4.1-8           scales_0.2.4          splines_3.1.2        
## [64] stats4_3.1.2          stringr_0.6.2         survival_2.37-7      
## [67] tcltk_3.1.2           tools_3.1.2           WGCNA_1.43           
## [70] XML_3.98-1.1          yaml_2.1.13
```

