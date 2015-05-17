## Heatmap examples



### Visualizing phyloseq objects

Also see the [phyloseq heatmap tutorial](http://joey711.github.io/phyloseq/plot_heatmap-examples.html).


```r
library(microbiome)
x <- download_microbiome("peerj32")$physeq
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

```r
plot_heatmap(x, sample.label = "sample")
```

```
## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
```



```r
library(microbiome)
x <- download_microbiome("dietswap")
```

```
## Downloading data set from O'Keefe et al. Nat. Comm. 6:6342, 2015 from Data Dryad: http://datadryad.org/resource/doi:10.5061/dryad.1mn1n
```

```r
bacteroidetes <- levelmap(NULL, "Phylum", "Genus", tax_table(x))$Bacteroidetes
x <- prune_taxa(bacteroidetes, x)
x <- subset_samples(x, group == "DI")
plot_heatmap(x, "NMDS", "bray", "sample", "Phylum")
```

```
## Error in .Method(..., na.last = na.last, decreasing = decreasing): argument 1 is not a vector
```




### Matrix visualization

Matrix heatmap


```r
library(microbiome)
data.peerj32 <- download_microbiome("peerj32")$physeq
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

```r
# Pick OTU Log10 matrix
x <- log10(otu_table(data.peerj32)@.Data)

# Visualize deviation of all bacteria from their population mean 
# (smaller: blue; higher: red):
# First Scale columns to zero mean and unit variance
# This puts the mean of each bacteria to zero 
# (less than mean is negative and blue; higher than mean is positive and red)
# having same variance for bacteria makes visual comparison easier too
x <- t(scale(t(x)))
tmp <- netresponse::plot_matrix(x, type = "twoway", mar = c(5, 14, 1, 1), cex.axis = 0.8)
```

![plot of chunk heatmap-matvisu-example](figure/heatmap-matvisu-example-1.png) 

Finding visually appealing order for rows and columns:


```r
hm <- heatmap(x) 
```

Then plot the same matrix with ordered rows and cols:


```r
tmp <- netresponse::plot_matrix(x[hm$rowInd, hm$colInd], type = "twoway", mar = c(5, 12, 1, 1), cex.axis = 0.5)
```

![plot of chunk heatmap-crosscorrelate3](figure/heatmap-crosscorrelate3-1.png) 


### Cross-correlating data sets

Cross-correlate columns of two data sets. This will return correlations, raw p-values, and fdr estimates (not strictly proper as the comparisons are not independent). Here we use robust biweight midcorrelation ('bicor') from the [WGCNA package](http://labs.genetics.ucla.edu/horvath/CoexpressionNetwork/Rpackages/WGCNA/). Keep only those elements that have at least only one significant correlation (n.signif):


```r
# Load example data
library(microbiome)
data.peerj32.otu <- download_microbiome("peerj32")$data$microbes
data.peerj32.lipids <- download_microbiome("peerj32")$data$lipids

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
head(correlation.table)
```

```
##                                   X1         X2 Correlation       p.adj
## 833      Ruminococcus gnavus et rel. TG.54.5..2   0.7207818 0.001738478
## 547      Ruminococcus gnavus et rel.   TG.52.5.   0.6996301 0.003192887
## 141 Eubacterium cylindroides et rel.   PC.40.3.  -0.6771286 0.003800575
## 144                     Helicobacter   PC.40.3.  -0.6838424 0.003800575
## 437      Ruminococcus gnavus et rel.   TG.50.4.   0.6852226 0.003800575
## 525      Ruminococcus gnavus et rel. TG.52.4..1   0.6716223 0.003800575
```

### Correlation heatmaps

Rearrange the data and plot the heatmap and mark significant correlations with stars to reproduce microbiota-lipidome heatmap from [this article](https://peerj.com/articles/32/): 


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



### Heatmaps

[Heatmaps](http://joey711.github.io/phyloseq/plot_heatmap-examples) and [Neatmaps](http://www.biomedcentral.com/1471-2105/11/45)


```r
library(microbiome)
data.dietswap <- download_microbiome("dietswap")

library(scales)
library(vegan)
plot_heatmap(physeq, taxa.label = "Phylum")
```

```
## Error in inherits(physeq, "formula"): object 'physeq' not found
```

```r
physeq.log <- physeq
```

```
## Error in eval(expr, envir, enclos): object 'physeq' not found
```

```r
trans <- as.matrix(t(scale(t(log10(1 + physeq.log@otu_table)))))
```

```
## Error in t(scale(t(log10(1 + physeq.log@otu_table)))): error in evaluating the argument 'x' in selecting a method for function 't': Error in t(log10(1 + physeq.log@otu_table)) : 
##   error in evaluating the argument 'x' in selecting a method for function 't': Error: object 'physeq.log' not found
```

```r
otu_table(physeq.log) <- otu_table(trans, taxa_are_rows=taxa_are_rows(physeq))
```

```
## Error in otu_table(trans, taxa_are_rows = taxa_are_rows(physeq)): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'trans' not found
```

```r
plot_heatmap(physeq.log, method = "NMDS", distance = "jaccard", low="blue", high="red")
```

```
## Error in inherits(physeq, "formula"): object 'physeq.log' not found
```

```r
#log transform sample counts
physeq.log <- transform_sample_counts(physeq, function(x) log10(x+1))
```

```
## Error in taxa_are_rows(physeq): error in evaluating the argument 'physeq' in selecting a method for function 'taxa_are_rows': Error: object 'physeq' not found
```

```r
#calculate ordination of log transformed sample counts
physeq.log.ord <- ordinate(physeq.log, "NMDS")
```

```
## Error in inherits(physeq, "formula"): object 'physeq.log' not found
```

```r
#split plot with log trandormed sample counts
p <- plot_ordination(physeq.log, physeq.log.ord, type = "split")
```

```
## Error in inherits(physeq, "phyloseq"): object 'physeq.log' not found
```

```r
p
```

![plot of chunk heatmap](figure/heatmap-1.png) 

