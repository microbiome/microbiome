<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{HITChip Atlas examples}
  %\usepackage[utf8]{inputenc}
-->


HITChip Atlas examples
======================

This tutorial provides step-by-step examples on exploratory analysis of
large-scale population-level microbiota profiling data.

### HITChip Atlas data set

An example data set from [Lahti et al. Nat. Comm. 5:4344,
2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)
contains large-scale profiling of 130 genus-like taxa across 1006 normal
western subjects. This data set is readily available for download from
the open [Data Dryad](http://doi.org/10.5061/dryad.pk75d) repository.

Load the HITChip Atlas microbiome profiling data in R.

    # Load Dryad tools
    library("rdryad") # Use the install.packages("rdryad") if package not available

    # Define the data URL
    url <- download_url('10255/dryad.64665')

    # Download the data
    data <- read.table(url, sep = "\t", row.names = 1, header = TRUE)

    # Fix some broken names from the original release..
    # ie. replace 'Clostridium..sensu.stricto.les' with 'Clostridiales'
    colnames(data) <- gsub("Clostridium..sensu.stricto.les", "Clostridiales", colnames(data))

Load the HITChip Atlas metadata in R. Note that some individuals have
multiple time points.

    url <- download_url('10255/dryad.64666')
    meta <- read.table(url, sep = "\t", row.names = 1, header = TRUE)

    # Add SampleIDs as a separate column from rownames
    meta$SampleID <- rownames(meta)

    # Order BMI groups in correct order
    # (see README at http://datadryad.org/resource/doi:10.5061/dryad.pk75d for details)
    meta$BMI_group <- factor(meta$BMI_group, levels = c("underweight", "lean", "overweight", "obese", "severeobese", "morbidobese"))
    meta$SubjectID <- factor(meta$SubjectID)

### Abundance histograms

Different sample sets have different population distributions in
microbial abundance. It is also important to consider whether to use
absolute or logarithmic abundances!

    # Load microbiome package
    library(microbiome)  

    # Load tools
    library(dplyr)

    # 1. List all samples (all time points and DNA extraction methods)
    all.samples <- meta$SampleID

    # 2. List samples at time point 0 that have specific DNA extraction method 
    rbb.samples <- filter(meta, Time == "0" & DNA_extraction_method == "r")$SampleID

    # Visualize
    #tax <- "Prevotella.melaninogenica.et.rel."
    tax <- "Bifidobacterium"
    d <- data[all.samples, tax]
    par(mfrow = c(1, 2))
    plot(density(d), main = paste(tax, "(All samples)"), xlab = "Abundance (Absolute HITChip signal)")
    plot(density(log10(d)), main = paste(tax, "(All samples)"), xlab = "Abundance (Log10 HITChip signal)")

![](Atlas_files/figure-markdown_strict/hist-1.png)

    d <- data[rbb.samples, tax]
    par(mfrow = c(1, 2))
    plot(density(d), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Absolute HITChip signal)")
    plot(density(log10(d)), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Log10 HITChip signal)")

![](Atlas_files/figure-markdown_strict/hist-2.png)

### Microbiota diversity

Diversity takes into account species richness and evenness ie. how
species abundances are distributed. We use here Shannon diversity.

    # Diversity using the vegan package
    # NOTE: data needs to be in absolute scale, not logarithmic
    library(vegan)
    di <- vegan::diversity(data, index = "shannon")

    # Diversity histogram across all samples
    hist(di, main = "Diversity histogram", xlab = "(Shannon) Diversity", ylab = "Population frequency")

![](Atlas_files/figure-markdown_strict/diversity-1.png)

Compare diversity with known background factors. By the way, we use here
all samples.

    par(mar = c(6, 4, 3, 1), mfrow = c(1, 2))
    boxplot(di ~ meta$BMI_group, las = 2, main = "Microbiota diversity vs. obesity")
    plot(meta$Age, di, main = "Microbiota diversity vs. Age", ylab = "Diversity", xlab = "Age (years)")

![](Atlas_files/figure-markdown_strict/diversitywithmetadata-1.png)

TASK: Try to use just single sample per subject (time point 0) and
perhaps a single DNA extraction method (r) - see above ?

### Plotting trends

Plot subject age versus phylotype abundance with smoothed confidence
intervals:

    library(microbiome)
    library(sorvi)
    df <- data.frame(Age = meta[rbb.samples, "Age"], Diversity = di[rbb.samples])
    p <- sorvi::regression_plot(Diversity~Age, df, shade = TRUE, mweight = TRUE, verbose = FALSE)
    print(p)

![](Atlas_files/figure-markdown_strict/visu-example3-1.png)

### Relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the input
data set needs to be in absolute scale (not logarithmic).

    rel <- relative.abundance(t(data))

    # Rearrange the data for ggplot visualization tools
    library(reshape)
    dfm <- melt(rel)
    colnames(dfm) <- c("Taxon", "SampleID", "RelativeAbundance")

    # Provide barplot visualizations of relative abundances for some randomly selected samples
    library(ggplot2)
    dfmf <- filter(dfm, SampleID %in% c("Sample-1", "Sample-2", "Sample-3", "Sample-4", "Sample-5"))
    p <- ggplot(dfmf, aes(x = SampleID, y = RelativeAbundance, fill = Taxon))
    p <- p + geom_bar(position = "stack", stat = "identity")
    print(p)

![](Atlas_files/figure-markdown_strict/diversity-example6-1.png)

    # Also note that taking relative abundances likely changes the abundance histograms
    theme_set(theme_bw(20))
    p <- ggplot(filter(dfm, Taxon == "Prevotella.melaninogenica.et.rel."), aes(x = 100*RelativeAbundance))
    p <- p + geom_density(fill = "darkgray")
    p <- p + scale_x_log10()
    p <- p + xlab("Relative Abundance (%)")
    print(p)

![](Atlas_files/figure-markdown_strict/diversity-example6-2.png)

### Principal component analysis (PCA)

Visualize deviation of all bacteria from their population mean (smaller:
blue; higher: red):

    # Let us focus on the most abundant and prevalent bacteria
    # that are seen in >1% (>0.01) relative abundance in 
    # >20% (>0.2) of the subjects
    prevalent.taxa <- names(which(prevalence(t(rel), 0.01, sort = TRUE) > 0.2))

    # Project data on 2D display with PCA (visualize subjects based on 20 random features)
    set.seed(423542)
    proj <- microbiome::project.data(log10(data[, prevalent.taxa]), type = "PCA")

    # Visualize
    p <- densityplot(proj, col = meta$DNA_extraction_method, legend = T)
    print(p)

![](Atlas_files/figure-markdown_strict/density-1.png)

    # Now do the same with RBB extracted samples only
    # Project data on 2D display with PCA (visualize subjects based on 20 random features)
    set.seed(4235423)
    proj <- microbiome::project.data(log10(data[rbb.samples, prevalent.taxa]), type = "PCA")

    # Visualize with DNA extraction method (now all samples have the same DNA extraction)
    p <- densityplot(proj, col = meta[rbb.samples, "DNA_extraction_method"], legend = T)
    print(p)

![](Atlas_files/figure-markdown_strict/density-2.png)

    # Visualize with low/high Prevotella
    # This shows that Prevotella (color) has ecosystem-level impact on microbiota composition
    #high.prevotella <- log10(data[rbb.samples, "Prevotella.melaninogenica.et.rel."]) > 4
    prevotella.abundance  <- log10(data[rbb.samples, "Prevotella.melaninogenica.et.rel."]) 
    p <- densityplot(proj, col = prevotella.abundance, legend = T)
    print(p)

![](Atlas_files/figure-markdown_strict/density-3.png)

PCA with ggplot2 - the above example gives a shortcut for the following:

    # Arrange projected data onto a data frame
    coms <- intersect(rownames(proj), rownames(meta))
    df <- as.data.frame(cbind(proj[coms,], meta[coms,]))
    names(df) <- c("x", "y", colnames(meta))

    # Construct the figure with ggplot2
    library(ggplot2)
    theme_set(theme_bw(15))
    p <- ggplot(df) 

    # Add densities
    p <- p + stat_density2d(aes(x = x, y = y, fill=..density..), geom="raster", stat_params = list(contour = F), geom_params = list()) 
    p <- p + scale_fill_gradient(low="white", high="black") 

    # Add points
    p <- p + geom_point(aes(x = x, y = y, color = Sex), size = 1.5) 

    # Add labels
    p <- p + xlab("PCA 1") + ylab("PCA 2") + ggtitle("Density plot")
    p <- p + scale_colour_discrete(name = "Sex")

    # Plot the figure
    print(p)

![](Atlas_files/figure-markdown_strict/density2-1.png)

### Licensing and Citations

This work can be freely used, modified and distributed under the
[Two-clause FreeBSD license](http://en.wikipedia.org/wiki/BSD_licenses).

Kindly cite the work as 'Leo Lahti and Gerben Hermes (2015). HITChip
Atlas tutorial. URL: <http://microbiome.github.com>'.

### Session info

This vignette was created with

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
    ## [1] parallel  stats     graphics  grDevices utils     datasets  methods  
    ## [8] base     
    ## 
    ## other attached packages:
    ##  [1] microbiome_0.99.34   ggplot2_1.0.0        sorvi_0.7.13        
    ##  [4] dplyr_0.3.0.2        rdryad_0.1.1         AnnotationDbi_1.26.1
    ##  [7] GenomeInfoDb_1.0.2   Biobase_2.24.0       BiocGenerics_0.10.0 
    ## [10] RSQLite_1.0.0        DBI_0.3.1            reshape_0.8.5       
    ## [13] vegan_2.2-1          lattice_0.20-29      permute_0.8-3       
    ## [16] e1071_1.6-4          devtools_1.7.0       knitcitations_1.0.5 
    ## [19] rmarkdown_0.3.10    
    ## 
    ## loaded via a namespace (and not attached):
    ##  [1] acepack_1.3-3.3       ape_3.1-4             assertthat_0.1       
    ##  [4] bibtex_0.4.0          class_7.3-11          cluster_1.15.3       
    ##  [7] codetools_0.2-9       colorspace_1.2-4      df2json_0.0.2        
    ## [10] digest_0.6.4          doParallel_1.0.8      dynamicTreeCut_1.62  
    ## [13] evaluate_0.5.5        fastcluster_1.1.15    foreach_1.4.2        
    ## [16] foreign_0.8-61        formatR_1.0           Formula_1.1-2        
    ## [19] gdata_2.13.3          GO.db_2.14.0          grid_3.1.2           
    ## [22] gtable_0.1.2          gtools_3.4.1          Hmisc_3.14-5         
    ## [25] htmltools_0.2.6       httr_0.5              igraph_0.7.1         
    ## [28] impute_1.38.1         IRanges_1.22.10       iterators_1.0.7      
    ## [31] knitr_1.8             labeling_0.3          latticeExtra_0.6-26  
    ## [34] lazyeval_0.1.9        lubridate_1.3.3       magrittr_1.0.1       
    ## [37] MASS_7.3-37           Matrix_1.1-4          matrixStats_0.10.3   
    ## [40] memoise_0.2.1         mgcv_1.8-3            mixOmics_5.0-3       
    ## [43] munsell_0.4.2         nlme_3.1-118          nnet_7.3-8           
    ## [46] OAIHarvester_0.1-7    pheatmap_0.7.7        plyr_1.8.1           
    ## [49] preprocessCore_1.26.1 proto_0.3-10          RColorBrewer_1.0-5   
    ## [52] Rcpp_0.11.3           RCurl_1.95-4.3        RefManageR_0.8.45    
    ## [55] reshape2_1.4.1        RGCCA_2.0             rgl_0.95.1158        
    ## [58] rjson_0.2.15          RJSONIO_1.3-0         R.methodsS3_1.6.1    
    ## [61] rpart_4.1-8           scales_0.2.4          splines_3.1.2        
    ## [64] stats4_3.1.2          stringr_0.6.2         survival_2.37-7      
    ## [67] tools_3.1.2           WGCNA_1.43            XML_3.98-1.1         
    ## [70] yaml_2.1.13
