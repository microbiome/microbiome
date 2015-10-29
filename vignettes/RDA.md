### RDA analysis and visualization. 

Load the package and example data:


```r
library(microbiome)
# Data from https://peerj.com/articles/32/
pseq <- download_microbiome("peerj32")$physeq
```

### Standard RDA 

Standard RDA for microbiota profiles versus the 'time' variable from 
sample metadata:


```r
rdatest <- rda_physeq(pseq, "time")
```

### RDA significance test


```r
permutest(rdatest) 
```

```
## Error in eval(expr, envir, enclos): could not find function "permutest"
```

### RDA visualization

Visualizing the standard RDA output:


```r
plot(rdatest, choices = c(1,2), type = "points", pch = 15, scaling = 3, cex = 0.7, col = meta$time)
points(rdatest, choices = c(1,2), pch = 15, scaling = 3, cex = 0.7, col = meta$time)
pl <- ordihull(rdatest, meta$time, scaling = 3, label = TRUE)
```

```
## Error in eval(expr, envir, enclos): could not find function "ordihull"
```

```r
title("RDA")
```

![plot of chunk rda4](figure/rda4-1.png) 

### Bagged RDA

Fitting bagged RDA on a phyloseq object:


```r
res <- bagged_rda(pseq, "group", sig.thresh=0.05, nboot=100)
```

```
## Error in eval(expr, envir, enclos): could not find function "bagged_rda"
```

Visualizing bagged RDA:


```r
plot_bagged_rda(res)
```

```
## Error in abs(Bag.res$loadings): non-numeric argument to mathematical function
```


### Controlling confounding variables with RDA

For more complex scenarios, use the vegan package directly:


```r
# Pick microbiota profiling data from the phyloseq object
otu <- otu_table(pseq)@.Data

# Sample annotations
meta <- sample_data(pseq)

# RDA with confounders
rdatest2 <- rda(t(otu) ~ meta$time + Condition(meta$subject + meta$gender))
```

```
## Error in eval(expr, envir, enclos): could not find function "rda"
```



