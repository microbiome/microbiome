### Oligoprofile heatmap

Reproduce and modify the oligo heatmap of the
[run.profiling.script](profiling). Plots heatmap of the oligo profiles
together with phylotype groups and sample clusters. Reload and plot
preprocessed data from the [run.profiling.script](profiling) by using
the [read.profiling](reading) function, assuming you have stored the
output in directory "datadirectory/". Start by reading the data:

    # Define data directory (here: simulated data directory)
    data.directory <- system.file("extdata", package = "microbiome")

    # Read Oligo level data in original domain
    oligo.matrix.log10.simulated <- read.profiling(level = "oligo", 
                              data.dir = data.directory, log10 = TRUE)

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/oligoprofile.tab
    ## Logarithmizing the data

    # Read Oligo-phylogeny mapping table (two methods):
    phylogeny.info <- GetPhylogeny("HITChip", "full")

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.full.tab

Save the image in a file (too large to be opened directly in R). To
prevent ordering of the rows, use hclust.method = NULL in the function
call:

    library(microbiome)
    ppcm <- 150
    png(filename = "oligoprofileClustering.png", 
         width = max(trunc(ppcm*21), 
                     trunc(ppcm*21*ncol(oligo.matrix.log10.simulated)/70)), 
         height = trunc(ppcm*29.7))


    library(HITChipDB)
    tmp <- PlotPhylochipHeatmap(oligo.matrix.log10.simulated, 
                    phylogeny.info, level = "L1", metric = "pearson")
    dev.off()
