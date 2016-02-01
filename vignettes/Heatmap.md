## Heatmaps

### phyloseq heatmaps

Download example data and plot heatmap with phyloseq tools:


```r
library(microbiome)
data("dietswap")
pseq <- dietswap

library(phyloseq)
plot_heatmap(pseq, sample.label = "sample")
```

![plot of chunk heatmap-phyloseq1](figure/heatmap-phyloseq1-1.png)

Pick subset of the data and plot ordered heatmap:


```r
bacteroidetes <- levelmap(from = "Phylum", to = "Genus",
	      	 	   tax.table = tax_table(pseq))$Bacteroidetes
pseq2 <- prune_taxa(bacteroidetes, pseq)
pseq2 <- subset_samples(pseq2, group == "DI")

# Plot heatmap ordered with NMDS
plot_heatmap(pseq2, method = "NMDS", distance = "bray",
	     sample.label = "sample", taxa.label = "Genus")
```

![plot of chunk heatmap-phyloseq2](figure/heatmap-phyloseq2-1.png)


### Matrix heatmaps

Alternatively, pick abundance matrix separately and use matrix
visualization tools. Z-transforming OTUs ie. visualize deviation of
all bacteria from their population mean (smaller: blue; higher: red):


```r
# Z transform
pseqz <- transform_phyloseq(pseq2, "Z", "OTU")

# Pick OTU table
x <- otu_table(pseqz)@.Data

# Plot heatmap
tmp <- netresponse::plot_matrix(x, type = "twoway", mar = c(5, 14, 1, 1))
```

![plot of chunk heatmap-matvisu-example](figure/heatmap-matvisu-example-1.png)


Find visually appealing order for rows and columns with Neatmap sorting:


```r
# Use the original phyloseq object for ordering as negative values are not allowed
data(peerj32)
x <- peerj32$microbes # Matrix
xz <- scale(x) # Z-transformed matrix

# Sorting rows (or columns) if just the sample order is needed rather than the matrix
sorted.rows <- neatsort(xz, "rows", method = "NMDS", distance = "euclidean") 

# Or just arrange the matrix directly
xo <- neat(xz, method = "NMDS", distance = "euclidean") # Sorted matrix
tmp <- netresponse::plot_matrix(t(xo), type = "twoway", mar = c(5, 12, 1, 1))
```

![plot of chunk neatmap3](figure/neatmap3-1.png)



### Cross-correlating data sets

Cross-correlate columns of two data sets from related to microbiome and blood serum lipids associations ([PeerJ 1:e32](https://peerj.com/articles/32/)).

The function returns correlations, raw p-values, and fdr estimates (not strictly proper as the comparisons are not independent). Here robust biweight midcorrelation ('bicor') from the [WGCNA package](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/). Keep only those elements that have at least only one significant correlation (n.signif):


```r
# Load example data 
library(microbiome)
data(peerj32)
data.peerj32.otu <- peerj32$microbes 
data.peerj32.lipids <- peerj32$lipids 

# Define data sets to cross-correlate
# OTU Log10 matrix # Microbiota (44 samples x 130 bacteria)
x <- log10(data.peerj32.otu)
y <- as.matrix(data.peerj32.lipids) # Lipids (44 samples x 389 lipids)

# Cross correlate data sets
correlations <- cross.correlate(x, y, method = "bicor", mode = "matrix", p.adj.threshold = 0.05, n.signif = 1)
```

Arrange the results in handy table format: 


```r
correlation.table <- cmat2table(correlations)
kable(head(correlation.table))
```



|    |X1                               |X2         | Correlation|     p.adj|
|:---|:--------------------------------|:----------|-----------:|---------:|
|833 |Ruminococcus gnavus et rel.      |TG(54:5).2 |   0.7207818| 0.0017385|
|547 |Ruminococcus gnavus et rel.      |TG(52:5)   |   0.6996301| 0.0031929|
|141 |Eubacterium cylindroides et rel. |PC(40:3)   |  -0.6771286| 0.0038006|
|144 |Helicobacter                     |PC(40:3)   |  -0.6838424| 0.0038006|
|437 |Ruminococcus gnavus et rel.      |TG(50:4)   |   0.6852226| 0.0038006|
|525 |Ruminococcus gnavus et rel.      |TG(52:4).1 |   0.6716223| 0.0038006|

### Correlation heatmaps

Rearrange the data and plot the heatmap and mark significant correlations with stars to reproduce microbiota-lipidome heatmap from [this article](https://peerj.com/articles/32/) (the ordering of rows and columns may be different): 


```r
p <- correlation.heatmap(correlation.table, "X1", "X2", fill = "Correlation", star = "p.adj", p.adj.threshold = 0.05) 
```

```r
print(p)
```

![plot of chunk heatmap-example-stars3](figure/heatmap-example-stars3-1.png)


### Heatmaps with ggplot2

The above examples provide handy shortcuts for heatmap visualization. You can also directly modify the ggplot2 routines. This time, let us set q-value threshold also for cell coloring: 


```r
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
```

![plot of chunk heatmap-example-stars](figure/heatmap-example-stars-1.png)

### Heatmap with text

For detailed information, might be handy to print the actual values on
top of the heatmap:


```r
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
```

![plot of chunk heatmap-example-text](figure/heatmap-example-text-1.png)

### Links

For further examples, see [phyloseq tutorial](http://joey711.github.io/phyloseq/plot_heatmap-examples.html) and [Neatmaps](http://www.biomedcentral.com/1471-2105/11/45)
