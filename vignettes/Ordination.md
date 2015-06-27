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

Project data on the first two PCA axes:


```r
set.seed(423542)
x <- t(otu_table(pseq2)@.Data)
```

```
## Error in t(otu_table(pseq2)@.Data): error in evaluating the argument 'x' in selecting a method for function 't': Error in otu_table(pseq2) : 
##   error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq2' not found
```

```r
proj <- project.data(log10(x), type = "PCA") # MDS.classical etc.
```

```
## Error in nrow(amat): object 'x' not found
```

Visualize and highlight:


```r
# Highlight gender
p <- densityplot(proj, col = sample_data(pseq2)$gender, legend = T)
```

```
## Error in as.vector(x, mode): cannot coerce type 'closure' to vector of type 'any'
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```

```r
# Highlight low/high Prevotella subjects
prevotella.abundance  <- log10(x[, "Prevotella melaninogenica et rel."]) 
```

```
## Error in eval(expr, envir, enclos): object 'x' not found
```

```r
p <- densityplot(proj, col = prevotella.abundance, legend = T)
```

```
## Error in as.vector(x, mode): cannot coerce type 'closure' to vector of type 'any'
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```

Projection with sample names:


```r
ggplot(aes(x = Comp.1, y = Comp.2, label = rownames(proj)), data = proj) + geom_text(size = 2)
```

```
## Error in eval(expr, envir, enclos): could not find function "ggplot"
```


### Unifrac with PCoA

See [phyloseq tutorial](http://joey711.github.io/phyloseq/plot_ordination-examples.html). 


### NMDS with Bray-Curtis distances


```r
nmds <- ordinate(pseq, "NMDS", "bray")
```

```
## Run 0 stress 0.09247795 
## Run 1 stress 0.1530495 
## Run 2 stress 0.143999 
## Run 3 stress 0.1531471 
## Run 4 stress 0.1419316 
## Run 5 stress 0.1179732 
## Run 6 stress 0.1495741 
## Run 7 stress 0.1453781 
## Run 8 stress 0.1529878 
## Run 9 stress 0.1410185 
## Run 10 stress 0.1415313 
## Run 11 stress 0.1311486 
## Run 12 stress 0.1508441 
## Run 13 stress 0.1217562 
## Run 14 stress 0.1477764 
## Run 15 stress 0.1526569 
## Run 16 stress 0.1446284 
## Run 17 stress 0.1492881 
## Run 18 stress 0.1553502 
## Run 19 stress 0.1527892 
## Run 20 stress 0.1503915
```

Ordinate the taxa in NMDS plot with Bray-Curtis distances


```r
pseq.ord <- ordinate(pseq2, "NMDS", "bray")
```

```
## Error in inherits(physeq, "formula"): object 'pseq2' not found
```

```r
p <- plot_ordination(pseq2, pseq.ord, type = "taxa", color = "Genus", title = "taxa")
```

```
## Error in inherits(physeq, "phyloseq"): object 'pseq2' not found
```

```r
print(p)
```

```
## Error in print(p): object 'p' not found
```

Grouping the plots by Phylum


```r
p + facet_wrap(~Phylum, 5)
```

```
## Error in eval(expr, envir, enclos): object 'p' not found
```



### Multidimensional scaling (MDS)


```r
plot_ordination(pseq, ordinate(pseq, "MDS"), color = "group") + geom_point(size = 5)
```

```
## Error in eval(expr, envir, enclos): could not find function "geom_point"
```


### Canonical correspondence analysis (CCA)

With samples:


```r
p <- plot_ordination(pseq, ordinate(pseq, "CCA"), type = "samples", color = "gender")
p + geom_point(size = 5) + geom_polygon(aes(fill = gender))
```

```
## Error in eval(expr, envir, enclos): could not find function "geom_point"
```

With taxa:


```r
plot_ordination(pseq, ordinate(pseq, "CCA"), type = "taxa", color = "Phylum") + 
    geom_point(size = 4)
```

```
## Error in eval(expr, envir, enclos): could not find function "geom_point"
```


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






