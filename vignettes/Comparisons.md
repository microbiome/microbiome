
## Group-wise comparisons

Read example data


```r
# Example data
library(microbiome)
pseq <- download_microbiome("dietswap")
```

### Comparing of two or more groups with a parametric test (linear model; ANOVA)


```r
# Define sample groups (gender + treatment group)
meta <- sample_data(pseq)
```

```
## Error in eval(expr, envir, enclos): could not find function "sample_data"
```

```r
group <- meta$group
```

```
## Error in eval(expr, envir, enclos): object 'meta' not found
```

```r
names(group) <- meta$sample
```

```
## Error in eval(expr, envir, enclos): object 'meta' not found
```

```r
# Calculate 1-way ANOVA p-values for the multi-group comparison
pvals <- check_anova(pseq, group, p.adjust.method = "BH", sort = TRUE)
```

```
## Error in eval(expr, envir, enclos): object 'group' not found
```


### Wilcoxon test (two-group comparisons)

If the data remarkably violates Gaussian assumptions, you need to use
non-parametric test. Wilcoxon is one option for two group
comparison. Here we compare males and females in the example data. 


```r
pval <- check_wilcoxon(pseq, "gender")
```

```
## Warning in check_foldchange(x, group, paired = paired): Converting the
## grouping variable gender into a factor.
```

```
## Error in check_foldchange(x, group, paired = paired): check_foldchange is currently implemented only for two-group comparisons. The selected variable gender has 0 levels:
```


