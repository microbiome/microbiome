## Ordination examples

Some examples with HITChip data. See also [phyloseq ordination tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html).

Load example data:


```r
library(microbiome)
library(phyloseq)
library(ggplot2)
data("atlas1006")
pseq <- atlas1006

# Convert signal to relative abundances
pseq.rel <- transform_phyloseq(pseq, "relative.abundance")

# Pick OTUs that are present with >1 percent relative abundance 
# in >10 percent of the samples
pseq2 <- filter_taxa(pseq.rel, function(x) sum(x > 1) > (0.1*nsamples(pseq.rel)), TRUE)
```


### Sample ordination

Project the samples with the given method and distance:


```r
set.seed(423542)
library(phyloseq)
pseq.ord <- ordinate(pseq2, "NMDS", "bray")
# Just pick the projected data (first two columns + metadata)
proj <- plot_ordination(pseq2, pseq.ord, justDF = T)
```


Visualize and highlight. In addition to densityplot, see plot_ordn from the microbiome package and plot_ordination from the phyloseq package.


```r
# Highlight gender
library(microbiome)
p <- densityplot(proj[, 1:2], col = proj$gender, legend = T)
print(p)
```

![plot of chunk ordination4](figure/ordination4-1.png)

```r
# Highlight low/high Prevotella subjects
prevotella.abundance  <- as.vector(log10(otu_table(pseq2)["Prevotella melaninogenica et rel.",]) )
p <- densityplot(proj[, 1:2], col = prevotella.abundance, legend = T)
print(p)
```

![plot of chunk ordination4](figure/ordination4-2.png)

Projection with sample names:


```r
ax1 <- names(proj)[[1]]
ax2 <- names(proj)[[2]]
ggplot(aes_string(x = ax1, y = ax2, label = "sample"), data = proj) + geom_text(size = 2)
```

![plot of chunk visu-example2](figure/visu-example2-1.png)


Ordinate the taxa in NMDS plot with Bray-Curtis distances


```r
p <- plot_ordination(pseq2, pseq.ord, type = "taxa", color = "Phylum", title = "Taxa ordination")
print(p)
```

![plot of chunk pca-ordination21](figure/pca-ordination21-1.png)

Grouping the plots by Phylum


```r
p + facet_wrap(~Phylum, 5)
```

![plot of chunk pca-ordination22](figure/pca-ordination22-1.png)


### Multidimensional scaling (MDS / PCoA)


```r
plot_ordination(pseq, ordinate(pseq, "MDS"), color = "DNA_extraction_method") + geom_point(size = 5)
```

![plot of chunk ordinate23](figure/ordinate23-1.png)

### RDA

See a separate page on [RDA](RDA.md).


### Canonical correspondence analysis (CCA)

With samples:


```r
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "samples", color = "gender")
p + geom_point(size = 5)
```

![plot of chunk ordinate24a](figure/ordinate24a-1.png)

With taxa:


```r
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "taxa", color = "Phylum")
p <- p + geom_point(size = 4)
print(p)
```

![plot of chunk ordinate24b](figure/ordinate24b-1.png)


### Split plot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "split", shape = "gender", 
    color = "Phylum", label = "gender")
```

![plot of chunk ordinate25](figure/ordinate25-1.png)


### Ordination biplot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "biplot", color = "Phylum")
```

![plot of chunk ordinate26](figure/ordinate26-1.png)






