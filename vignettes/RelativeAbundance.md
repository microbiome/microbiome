### Relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the input
data set needs to be in absolute scale (not logarithmic).

    library(microbiome)
    # Define data path (you can replace data.directory with your own path)
    data.directory <- system.file("extdata", package = "microbiome")
    x <- read.profiling(level = "L1", method = "frpa", 
                          data.dir = data.directory, 
                          log10 = FALSE)  

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/L1-frpa.tab

    rel <- relative.abundance(x)

    # Rearrange the data for ggplot visualization tools
    library(reshape)
    dfm <- melt(rel)
    colnames(dfm) <- c("Taxon", "SampleID", "RelativeAbundance")

    # Provide barplot visualizations of relative abundances for some randomly selected samples
    library(ggplot2)
    suppressMessages(library(dplyr))
    dfmf <- filter(dfm, SampleID %in% c("Sample.1", "Sample.2", "Sample.3", "Sample.4", "Sample.5"))
    p <- ggplot(dfmf, aes(x = SampleID, y = RelativeAbundance, fill = Taxon))
    p <- p + geom_bar(position = "stack", stat = "identity")
    print(p)

![](figure/diversity-example6-1.png)

    # Also note that taking relative abundances likely changes the abundance histograms
    tax <- "Bacilli"
    theme_set(theme_bw(20))
    p <- ggplot(filter(dfm, Taxon == tax), aes(x = 100*RelativeAbundance))
    p <- p + geom_density(fill = "darkgray")
    p <- p + scale_x_log10()
    p <- p + xlab("Relative Abundance (%)")
    p <- p + ggtitle(tax)
    print(p)

![](figure/diversity-example6-2.png)
