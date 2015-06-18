## Ordination examples

Here some examples with HITChip data. See also [phyloseq ordination tutorial](joey711.github.io/phyloseq/plot_ordination-examples.html).


Load example data:


```r
library(microbiome)
pseq.atlas <- download_microbiome("atlas1006")

# Define the data set to be used in the examples:
pseq <- pseq.atlas

# Convert signal to relative abundances
pseq.rel <- transform_sample_counts(pseq, function (x) {x/sum(x)})

# Pick OTUs that are present with >1 percent relative abundance 
# in >10 percent of the samples
pseq2 <- filter_taxa(pseq.rel, function(x) sum(x > .01) > (0.1*nsamples(pseq.rel)), TRUE)
```


### Principal component analysis (PCA)


```r
# Project data on 2D display with PCA
set.seed(423542)
x <- t(otu_table(pseq2)@.Data)
proj <- project.data(log10(x), type = "PCA") # MDS.classical etc.
```

Visualize and highlight:


```r
# Highlight gender
p <- densityplot(proj, col = sample_data(pseq2)$gender, legend = T)
print(p)
```

![plot of chunk ordination4](figure/ordination4-1.png) 

```r
# Highlight low/high Prevotella subjects
prevotella.abundance  <- log10(x[, "Prevotella melaninogenica et rel."]) 
p <- densityplot(proj, col = prevotella.abundance, legend = T)
print(p)
```

![plot of chunk ordination4](figure/ordination4-2.png) 

Projection with sample names:


```r
ggplot(aes(x = Comp.1, y = Comp.2, label = rownames(proj)), data = proj) + geom_text(size = 2)
```

![plot of chunk visu-example2](figure/visu-example2-1.png) 


### Unifrac with PCoA

See [phyloseq tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html). 


### NMDS with Bray-Curtis distances


```r
nmds <- ordinate(pseq, "NMDS", "bray")
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.2044329 
## Run 1 stress 0.2047176 
## ... procrustes: rmse 0.01353912  max resid 0.1388324 
## Run 2 stress 0.203759 
## ... New best solution
## ... procrustes: rmse 0.009482057  max resid 0.1384349 
## Run 3 stress 0.2053158 
## Run 4 stress 0.2168802 
## Run 5 stress 0.2049942 
## Run 6 stress 0.2063315 
## Run 7 stress 0.2089041 
## Run 8 stress 0.2087436 
## Run 9 stress 0.2049411 
## Run 10 stress 0.206135 
## Run 11 stress 0.2083954 
## Run 12 stress 0.2136372 
## Run 13 stress 0.2086734 
## Run 14 stress 0.2098592 
## Run 15 stress 0.2121382 
## Run 16 stress 0.2074855 
## Run 17 stress 0.2050563 
## Run 18 stress 0.2084222 
## Run 19 stress 0.2084397 
## Run 20 stress 0.2158551
```

Ordinate the taxa in NMDS plot with Bray-Curtis distances


```r
pseq.ord <- ordinate(pseq2, "NMDS", "bray")
```

```
## Square root transformation
## Wisconsin double standardization
## Run 0 stress 0.2014154 
## Run 1 stress 0.4170224 
## Run 2 stress 0.2175427 
## Run 3 stress 0.2120913 
## Run 4 stress 0.2149669 
## Run 5 stress 0.2205538 
## Run 6 stress 0.23016 
## Run 7 stress 0.2120833 
## Run 8 stress 0.2174832 
## Run 9 stress 0.2114263 
## Run 10 stress 0.2116162 
## Run 11 stress 0.2013954 
## ... New best solution
## ... procrustes: rmse 0.008975218  max resid 0.09184326 
## Run 12 stress 0.2081587 
## Run 13 stress 0.2284501 
## Run 14 stress 0.2142035 
## Run 15 stress 0.2039289 
## Run 16 stress 0.2141113 
## Run 17 stress 0.2041206 
## Run 18 stress 0.2362722 
## Run 19 stress 0.2102958 
## Run 20 stress 0.2085741
```

```r
p <- plot_ordination(pseq2, pseq.ord, type = "taxa", color = "Genus", title = "taxa")
print(p)
```

![plot of chunk pca-ordination21](figure/pca-ordination21-1.png) 

Grouping the plots by Phylum


```r
p + facet_wrap(~Phylum, 5)
```

![plot of chunk pca-ordination22](figure/pca-ordination22-1.png) 



### Multidimensional scaling (MDS)


```r
plot_ordination(pseq, ordinate(pseq, "MDS"), color = "group") + geom_point(size = 5)
```

![plot of chunk ordinate23](figure/ordinate23-1.png) 


### Canonical correspondence analysis (CCA)

With samples:


```r
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "samples", color = "gender")
p + geom_point(size = 5) + geom_polygon(aes(fill = gender))
```

![plot of chunk ordinate24a](figure/ordinate24a-1.png) 

With taxa:


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "taxa", color = "Phylum") + 
    geom_point(size = 4)
```

![plot of chunk ordinate24b](figure/ordinate24b-1.png) 


### Split plot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "split", color = "gender", 
    shape = "Phylum", label = "gender")
```

![plot of chunk ordinate25](figure/ordinate25-1.png) 


### Ordination biplot


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "biplot", shape = "Phylum")
```

![plot of chunk ordinate26](figure/ordinate26-1.png) 

See [phyloseq tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html). 




