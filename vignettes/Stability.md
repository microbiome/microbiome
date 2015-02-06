### Stability analysis within group of samples

Calculate stability as the average correlation between samples and their
mean for a given phylotypes vs. samples matrix:

    # Example data
    library(microbiome)
    data(peerj32)
    x <- peerj32$microbes

    # Estimate Stability
    stability <- estimate.stability(t(x))$stability

### Stability analysis in time

Calculate correlations between time points 1 and 2 for each subject,
then calculate stability as the average correlation of these
subject-specific correlations. Done as above, but instead of one matrix,
give two matrices as input. These are phylotypes vs. samples matrices
and correspond to time points 1 and 2, respectively. The rows
(phylotypes) and columns (subjects) should be in the same order in both
matrices.

### Version information

    sessionInfo()

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
    ##  [1] tcltk     grid      parallel  stats     graphics  grDevices utils    
    ##  [8] datasets  methods   base     
    ## 
    ## other attached packages:
    ##  [1] googleVis_0.5.8      HITChipDB_0.5.13     RPA_1.20.01         
    ##  [4] affydata_1.12.0      affy_1.42.3          RMySQL_0.10.1       
    ##  [7] ade4_1.6-2           limma_3.20.9         gdata_2.13.3        
    ## [10] netresponse_1.17.13  mclust_4.4           minet_3.20.2        
    ## [13] infotheo_1.2.0       Rgraphviz_2.8.1      graph_1.42.0        
    ## [16] ggplot2_1.0.0        knitr_1.9            dplyr_0.4.1         
    ## [19] sorvi_0.7.13         microbiome_0.99.35   AnnotationDbi_1.26.1
    ## [22] GenomeInfoDb_1.0.2   Biobase_2.24.0       BiocGenerics_0.10.0 
    ## [25] RSQLite_1.0.0        DBI_0.3.1            reshape_0.8.5       
    ## [28] vegan_2.2-1          lattice_0.20-29      permute_0.8-3       
    ## [31] e1071_1.6-4          rdryad_0.1.1         knitcitations_1.0.5 
    ## [34] rmarkdown_0.5.1     
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] acepack_1.3-3.3       affyio_1.32.0         ape_3.2              
    ##  [4] assertthat_0.1        bibtex_0.4.0          BiocInstaller_1.14.3 
    ##  [7] bitops_1.0-6          class_7.3-11          cluster_1.15.3       
    ## [10] codetools_0.2-9       colorspace_1.2-4      df2json_0.0.2        
    ## [13] digest_0.6.8          dmt_0.8.20            doParallel_1.0.8     
    ## [16] dynamicTreeCut_1.62   evaluate_0.5.5        fastcluster_1.1.15   
    ## [19] foreach_1.4.2         foreign_0.8-61        formatR_1.0          
    ## [22] Formula_1.2-0         GO.db_2.14.0          gtable_0.1.2         
    ## [25] gtools_3.4.1          Hmisc_3.14-6          htmltools_0.2.6      
    ## [28] httr_0.6.1            igraph_0.7.1          impute_1.38.1        
    ## [31] IRanges_1.22.10       iterators_1.0.7       labeling_0.3         
    ## [34] latticeExtra_0.6-26   lazyeval_0.1.10       lubridate_1.3.3      
    ## [37] magrittr_1.5          MASS_7.3-37           Matrix_1.1-5         
    ## [40] matrixStats_0.13.1    memoise_0.2.1         mgcv_1.8-3           
    ## [43] mixOmics_5.0-3        munsell_0.4.2         mvtnorm_1.0-2        
    ## [46] nlme_3.1-119          nnet_7.3-8            OAIHarvester_0.1-7   
    ## [49] pheatmap_0.7.7        plyr_1.8.1            preprocessCore_1.26.1
    ## [52] proto_0.3-10          qvalue_1.38.0         RColorBrewer_1.1-2   
    ## [55] Rcpp_0.11.4           RCurl_1.95-4.5        RefManageR_0.8.45    
    ## [58] reshape2_1.4.1        RGCCA_2.0             rgl_0.95.1201        
    ## [61] rjson_0.2.15          RJSONIO_1.3-0         R.methodsS3_1.6.1    
    ## [64] rpart_4.1-8           scales_0.2.4          splines_3.1.2        
    ## [67] stats4_3.1.2          stringr_0.6.2         survival_2.37-7      
    ## [70] tools_3.1.2           WGCNA_1.43            XML_3.98-1.1         
    ## [73] yaml_2.1.13           zlibbioc_1.10.0
