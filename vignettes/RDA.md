<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - rda}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
RDA analysis and visualization.
-------------------------------

Load the package and example data:

    library(microbiome)
    data(peerj32) # Data from https://peerj.com/articles/32/
    pseq <- peerj32$phyloseq

### Standard RDA

Standard RDA for microbiota profiles versus the given (here 'time')
variable from sample metadata:

    pseq.trans <- transform(pseq, "hell") # Hellinger transform
    rda.result <- rda_physeq(pseq.trans, "time", scale = TRUE)

    # Proportion explained by the contraints
    summary(rda.result)$constr.chi/summary(rda.result)$tot.chi

    ## [1] 0.01540884

### RDA visualization

Visualizing the standard RDA output:

    library(phyloseq)
    meta <- sample_data(pseq.trans)
    plot(rda.result, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = meta$time)
    points(rda.result, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = meta$time)
    library(vegan)
    pl <- ordihull(rda.result, meta$time, scaling = 3, label = TRUE)
    title("RDA")

![](RDA_files/figure-markdown_strict/rda4-1.png)

See also the RDA method in phyloseq::ordinate, which is calculated
without the formula.

### RDA significance test

    library(vegan)
    permutest(rda.result) 

    ## 
    ## Permutation test for rda 
    ## 
    ## Permutation: free
    ## Number of permutations: 99
    ##  
    ## Call: rda(formula = otu ~ annot, scale = scale, na.action =
    ## na.action)
    ## Permutation test for all constrained eigenvalues
    ## Pseudo-F:     0.6572996 (with 1, 42 Degrees of Freedom)
    ## Significance:     0.89

### Bagged RDA

Fitting bagged (bootstrap aggregated) RDA on a phyloseq object:

    res <- bagged_rda(pseq.trans, "group", sig.thresh=0.05, nboot=100)

Visualizing bagged RDA:

    plot_bagged_rda(res)

![](RDA_files/figure-markdown_strict/rda6-1.png)

### RDA with confounding variables

For more complex RDA scenarios, use the vegan package:

    # Pick microbiota profiling data from the phyloseq object
    otu <- abundances(pseq.trans)

    # Sample annotations
    meta <- sample_data(pseq.trans)

    # RDA with confounders
    rda.result2 <- rda(t(otu) ~ meta$time + Condition(meta$subject + meta$gender))
