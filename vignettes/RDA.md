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
    # Data from https://peerj.com/articles/32/
    data("peerj32")
    pseq <- peerj32$phyloseq

### Standard RDA

Standard RDA for microbiota profiles versus the given (here 'time')
variable from sample metadata:

    pseq.trans <- transform_phyloseq(pseq, "hell") # Hellinger transform
    rda.result <- rda_physeq(pseq.trans, "time", scale = TRUE)

    # Proportion explained by the contraints
    summary(rda.result)$constr.chi/summary(rda.result)$tot.chi

### RDA visualization

Visualizing the standard RDA output:

    library(phyloseq)
    meta <- sample_data(pseq.trans)
    plot(rda.result, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = meta$time)
    points(rda.result, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = meta$time)
    library(vegan)
    pl <- ordihull(rda.result, meta$time, scaling = 3, label = TRUE)
    title("RDA")

See also the RDA method in phyloseq::ordinate, which is calculated
without the formula.

### RDA significance test

    library(vegan)
    permutest(rda.result) 

### Bagged RDA

Fitting bagged (bootstrap aggregated) RDA on a phyloseq object:

    res <- bagged_rda(pseq.trans, "group", sig.thresh=0.05, nboot=100)

Visualizing bagged RDA:

    plot_bagged_rda(res)

### RDA with confounding variables

For more complex RDA scenarios, use the vegan package:

    # Pick microbiota profiling data from the phyloseq object
    otu <- abundances(pseq.trans)

    # Sample annotations
    meta <- sample_data(pseq.trans)

    # RDA with confounders
    rda.result2 <- rda(t(otu) ~ meta$time + Condition(meta$subject + meta$gender))
