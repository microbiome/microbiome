<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - heatmap}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Heatmaps for microbiome analysis
--------------------------------

See [Composition](Composition.html) page for further microbiota
composition heatmaps, as well as the [phyloseq
tutorial](http://joey711.github.io/phyloseq/plot_heatmap-examples.html)
and [Neatmaps](http://www.biomedcentral.com/1471-2105/11/45). Moreover,
the [aheatmap](http://nmf.r-forge.r-project.org/aheatmap.html) function
of the NMF package provides further high quality heatmap plotting
capabilities with row and column annotation color bars, clustering trees
and other useful features that are often missing from standard heatmap
tools in R.

Load some example data:

    library(microbiome) # Load libraries
    library(phyloseq)
    data(peerj32)
    pseq <- peerj32$phyloseq    # Rename data

    # Pick data subset (DI samples from Phylum Bacteroidetes)
    pseq2 <- pseq %>%
             subset_taxa(Phylum == "Bacteroidetes") %>%
             subset_samples(group == "LGG")

    # Z transformed abundance data
    pseqz <- microbiome::transform(pseq2, "Z")

### Matrix heatmaps

Visualize the Z-transformed abundance matrix

    # Plot the abundances heatmap
    heat(melt(abundances(pseqz)), "Var1", "Var2", "value")

![](Heatmap_files/figure-markdown_strict/hn122-1.png)

Find visually appealing order for rows and columns with the Neatmap
approach:

    # Sort the matrix rows and cols directly
    xo <- neat(abundances(pseqz), method = "NMDS", distance = "euclidean") 

    # Heatmap visualization
    heat(melt(xo), "Var1", "Var2", "value")

![](Heatmap_files/figure-markdown_strict/neatmap3-1.png)

    # or use a shortcut to sorting rows (or columns) if just the order was needed 
    sorted.rows <- neatsort(abundances(pseqz), "rows", method = "NMDS", distance = "euclidean") 

### Cross-correlating data sets

Cross-correlate columns of two data sets from related to microbiome and
blood serum lipids associations ([PeerJ
1:e32](https://peerj.com/articles/32/)).

The function returns correlations, raw p-values, and fdr estimates (not
strictly proper as the comparisons are not independent). Keep only those
elements that have at least only one significant correlation (n.signif):

    # Load example data 
    otu <- peerj32$microbes 
    lipids <- peerj32$lipids 

    # Define data sets to cross-correlate
    x <- log10(otu) # OTU Log10 (44 samples x 130 genera)
    y <- as.matrix(lipids) # Lipids (44 samples x 389 lipids)

    # Cross correlate data sets
    correlations <- associate(x, y, method = "spearman", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)

    # Or, alternatively, the same output is also available in a handy table format
    correlation.table <- associate(x, y, method = "spearman", mode = "table", p.adj.threshold = 0.05, n.signif = 1)

    kable(head(correlation.table))

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">X1</th>
<th align="left">X2</th>
<th align="right">Correlation</th>
<th align="right">p.adj</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>1648</td>
<td align="left">Ruminococcus gnavus et rel.</td>
<td align="left">TG(54:5).2</td>
<td align="right">0.7164958</td>
<td align="right">0.0022842</td>
</tr>
<tr class="even">
<td>384</td>
<td align="left">Moraxellaceae</td>
<td align="left">PC(40:3e)</td>
<td align="right">-0.6941863</td>
<td align="right">0.0029225</td>
</tr>
<tr class="odd">
<td>1829</td>
<td align="left">Uncultured Bacteroidetes</td>
<td align="left">TG(56:2).1</td>
<td align="right">-0.6987375</td>
<td align="right">0.0029225</td>
</tr>
<tr class="even">
<td>349</td>
<td align="left">Lactobacillus plantarum et rel.</td>
<td align="left">PC(40:3)</td>
<td align="right">-0.6877976</td>
<td align="right">0.0031520</td>
</tr>
<tr class="odd">
<td>1198</td>
<td align="left">Ruminococcus gnavus et rel.</td>
<td align="left">TG(52:5)</td>
<td align="right">0.6806216</td>
<td align="right">0.0037518</td>
</tr>
<tr class="even">
<td>264</td>
<td align="left">Moraxellaceae</td>
<td align="left">PC(38:4).1</td>
<td align="right">-0.6700504</td>
<td align="right">0.0038414</td>
</tr>
</tbody>
</table>

### Association heatmaps

Rearrange the data and plot the heatmap and mark significant
correlations with stars to reproduce microbiota-lipidome heatmap from
[Lahti et al. PeerJ (2013)](https://peerj.com/articles/32/) (the
ordering of rows and columns may be different):

    p <- heat(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 

    print(p)

![](Heatmap_files/figure-markdown_strict/heatmap-example-stars3-1.png)

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

![](Heatmap_files/figure-markdown_strict/heatmap-example-stars-1.png)

### Heatmap with text

For detailed information, might be handy to print the actual values on
top of the heatmap:

    theme_set(theme_bw(20))
    df <- correlation.table
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

![](Heatmap_files/figure-markdown_strict/heatmap-example-text-1.png)

### ggcorr

An alternative way to visualize correlation matrices is provided by the
[ggcorr package](https://briatte.github.io/ggcorr/). Note: this toy
example does not consider the compositionality effect in microbial
abundance correlations. See the package site for more detailed examples
and many more options.

    library(GGally)
    ggcorr(x[, 1:10], method = c("pairwise", "spearman"), nbreaks = 20, hjust = 0.75)
    ggcorr(x[, 1:10], method = c("pairwise", "spearman"), nbreaks = 20, geom = "circle")
    ggcorr(x[, 1:10], method = c("pairwise", "spearman"), nbreaks = 20, label = TRUE, label_alpha = TRUE)
    ggcorr(data = NULL, cor_matrix = cor(x[, 1:10], use = "everything"), low = "steelblue", mid = "white", high = "darkred", midpoint = 0)

<img src="Heatmap_files/figure-markdown_strict/ggcorr1-1.png" width="400px" /><img src="Heatmap_files/figure-markdown_strict/ggcorr1-2.png" width="400px" /><img src="Heatmap_files/figure-markdown_strict/ggcorr1-3.png" width="400px" /><img src="Heatmap_files/figure-markdown_strict/ggcorr1-4.png" width="400px" />
