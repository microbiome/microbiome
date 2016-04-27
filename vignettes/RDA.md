## RDA analysis and visualization. 

Load the package and example data:


```r
library(microbiome)
# Data from https://peerj.com/articles/32/
data("peerj32")
pseq <- peerj32$phyloseq
```

### Standard RDA 

Standard RDA for microbiota profiles versus the given (here 'time')
variable from sample metadata:


```r
# If x has zeroes we can use log(1 + x) transformation
pseq.log10 <- transform_phyloseq(pseq, "log10")
rda.result <- rda_physeq(pseq.log10, "time")
```

```
## Error in qr.fitted(Q, Xbar): 'qr' and 'y' must have the same number of rows
```

### RDA visualization

Visualizing the standard RDA output:


```r
library(phyloseq)
meta <- sample_data(pseq.log10)
plot(rda.result, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = meta$time)
```

```
## Error in plot(rda.result, choices = c(1, 2), type = "points", pch = 15, : object 'rda.result' not found
```

```r
points(rda.result, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = meta$time)
```

```
## Error in points(rda.result, choices = c(1, 2), pch = 15, scaling = 3, : object 'rda.result' not found
```

```r
library(vegan)
pl <- ordihull(rda.result, meta$time, scaling = 3, label = TRUE)
```

```
## Error in scores(ord, display = display, ...): object 'rda.result' not found
```

```r
title("RDA")
```

```
## Error in title("RDA"): plot.new has not been called yet
```

See also the RDA method in phyloseq ordinate function, which is calculated without the formula.


### RDA significance test


```r
library(vegan)
permutest(rda.result) 
```

```
## Error in permutest(rda.result): object 'rda.result' not found
```

### Bagged RDA

Fitting bagged (bootstrap aggregated) RDA on a phyloseq object:


```r
res <- bagged_rda(pseq.log10, "group", sig.thresh=0.05, nboot=100)
```

Visualizing bagged RDA:


```r
plot_bagged_rda(res)
```

![plot of chunk rda6](figure/rda6-1.png)


### RDA with confounding variables 

For more complex RDA scenarios, use the vegan package directly:


```r
# Pick microbiota profiling data from the phyloseq object
otu <- otu_table(pseq.log10)@.Data

# Sample annotations
meta <- sample_data(pseq.log10)

# RDA with confounders
rda.result2 <- rda(t(otu) ~ meta$time + Condition(meta$subject + meta$gender))
```



