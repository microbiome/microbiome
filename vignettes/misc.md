### Bagged RDA

Calculate Bagged RDA and visualize the results:

    # Modify the group vector into a facor as required by Bagged RDA
    y <- factor(annot$time); names(y) <- rownames(annot)

    # Bagged RDA
    Bag.res <- Bagged.RDA.Feature.Selection(t(l2), y, sig.thresh=0.05, nboot=100)

    # Visualize
    PlotBaggedRDA(Bag.res, y)

### Oligo heatmap

This reproduces the oligo-level heatmap from profiling script. If there
are problems, try to tune ppcm, figureratio and fontsize (see
help(add.heatmap) for details)

    library(microbiome)

    # Load Phylogeny
    phylogeny.info <- GetPhylogeny("HITChip")

    # Load oligo-level data

    # Replace data.directory here with your own profiling script output data directory
    data.directory <- system.file("extdata", package = "microbiome")

    oligodata <- read.profiling(level = "oligo", log10 = FALSE, data.dir = data.directory)

    # Produce the plot and save it to the working directory
    library(HITChipDB)

    ##      Package    LibPath                                            
    ## [1,] "affydata" "/home/antagomir/R/x86_64-pc-linux-gnu-library/3.1"
    ##      Item       Title                        
    ## [1,] "Dilution" "AffyBatch instance Dilution"

    hc.params <- add.heatmap(log10(oligodata), output.dir = ".", phylogeny.info = phylogeny.info)
