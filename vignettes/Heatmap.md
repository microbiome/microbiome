### Matrix visualization

Matrix heatmap

    # Load example data
    library(microbiome)
    data(peerj32)
    x <- peerj32$microbes

    # Log10
    x <- log10(x)

    # Transponse the matrix 
    x <- t(x)

    # Visualize deviation of all bacteria from their population mean (smaller: blue; higher: red):
    # First Scale columns to zero mean and unit variance
    # This puts the mean of each bacteria to zero 
    # (less than mean is negative and blue; higher than mean is positive and red)
    # having same variance for bacteria makes visual comparison easier too
    x <- scale(x)
    tmp <- netresponse::plot_matrix(x, type = "twoway", mar = c(5, 14, 1, 1), cex.axis = 0.8)

![](figure/matvisu-example-1.png)

Finding visually appealing order for rows and columns:

    hm <- heatmap(x) 

Then plot the same matrix with ordered rows and cols:

    tmp <- netresponse::plot_matrix(x[hm$rowInd, hm$colInd], type = "twoway", mar = c(5, 12, 1, 1), cex.axis = 0.5)

![](figure/heatmap-crosscorrelate3-1.png)

### Cross-correlating data sets

Cross-correlate columns of two data sets. This will return correlations,
raw p-values, and fdr estimates (not strictly proper as the comparisons
are not independent). Here we use robust biweight midcorrelation
('bicor') from the [WGCNA
package](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/).
Keep only those elements that have at least only one significant
correlation (n.signif):

    # Load example data
    library(microbiome)
    data(peerj32)
    x <- peerj32$microbes

    # Define data sets to cross-correlate
    dat1 <- peerj32$lipids # Lipids (44 samples x 389 lipids)
    dat2 <- log10(peerj32$microbes) # Microbiota (44 samples x 130 bacteria)

    correlations <- microbiome::cross.correlate(dat1, dat2, 
                                method = "bicor", mode = "matrix", 
                            n.signif = 1, p.adj.threshold = 0.05, 
                            p.adj.method = "fdr")

Arrange the results in handy table format:

    correlation.table <- microbiome::cmat2table(correlations)
    head(correlation.table)

    ##              X1                               X2 Correlation       p.adj
    ## 1100 TG(54:5).2      Ruminococcus gnavus et rel.   0.7207818 0.001738478
    ## 1087   TG(52:5)      Ruminococcus gnavus et rel.   0.6996301 0.003192887
    ## 479    PC(40:3) Eubacterium cylindroides et rel.  -0.6771286 0.003800575
    ## 656    PC(40:3)                     Helicobacter  -0.6838424 0.003800575
    ## 1082   TG(50:4)      Ruminococcus gnavus et rel.   0.6852226 0.003800575
    ## 1086 TG(52:4).1      Ruminococcus gnavus et rel.   0.6716223 0.003800575

### Correlation heatmaps

Rearrange the data and plot the heatmap and mark significant
correlations with stars to reproduce microbiota-lipidome heatmap from
[this article](https://peerj.com/articles/32/):

    p <- microbiome::correlation.heatmap(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 

    print(p)

![](figure/heatmap-example-stars3-1.png)

### Heatmaps with ggplot2

The above examples provide handy shortcuts for heatmap visualization.
You can also directly modify the ggplot2 routines. This time, let us set
q-value threshold also for cell coloring:

    # Order the rows and columns with levels argument if needed:
    correlation.table$X1 <- factor(correlation.table$X1, levels = unique(as.character(correlation.table$X1)))
    correlation.table$X2 <- factor(correlation.table$X2, levels = unique(as.character(correlation.table$X2)))

    # Set black-and-white theme
    library(ggplot2)
    theme_set(theme_bw())

    # Pick only the correlations with q<0.05
    # Note: this will leave other cells empty
    library(dplyr)
    subtable <- filter(correlation.table, p.adj < 0.05)

    # Arrange the figure
    p <- ggplot(subtable, aes(x = X1, y = X2, fill = Correlation))
    p <- p + geom_tile() 
    p <- p + scale_fill_gradientn("Correlation", 
                           breaks = seq(from = -1, to = 1, by = 0.2), 
                       colours = c("darkblue", "blue", "white", "red", "darkred"), 
                       limits = c(-1,1)) 

    # Polish texts
    p <- p + theme(axis.text.x=element_text(angle = 90))
    p <- p + xlab("") + ylab("")

    # Mark the most significant cells with stars
    p <- p + geom_text(data = subset(correlation.table, p.adj < 0.02), 
                   aes(x = X1, y = X2, label = "+"), col = "white", size = 5)

    # Plot
    print(p)

![](figure/heatmap-example-stars-1.png)

### Heatmap with text

For detailed information, might be handy to print the actual values on
top of the heatmap:

    theme_set(theme_bw(20))
    df <- microbiome::cmat2table(correlations)
    df$X1 <- factor(df$X1)
    df$X2 <- factor(df$X2)
    p <- ggplot(df, aes(X1, X2, group=X2)) 
    p <- p + geom_tile(aes(fill = Correlation)) 
    p <- p + geom_text(aes(fill = df$Correlation, label = round(df$Correlation, 1)), size = 2) 
    p <- p + scale_fill_gradientn("Correlation", 
                          breaks = seq(from = -1, to = 1,  by = 0.25), 
                          colours = c("blue", "white", "red"), 
                      limits = c(-1, 1))
    p <- p + theme(axis.text.x = element_text(angle = 90, vjust = 0.5)) 
    p <- p + xlab("") + ylab("")
    print(p)

![](figure/heatmap-example-text-1.png)
