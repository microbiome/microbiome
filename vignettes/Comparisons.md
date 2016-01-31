
## Group-wise comparisons

Read example data


```r
# Example data
library(microbiome)
pseq <- download_microbiome("dietswap")
```

```
## Error in curl::curl_fetch_memory(url, handle = handle): Timeout was reached
```

### Comparing of two or more groups with a parametric test (linear model; ANOVA)


```r
# Define sample groups (gender + treatment group)
meta <- sample_data(pseq)
```

```
## Error in sample_data(pseq): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq' not found
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
## Error in check_anova(pseq, group, p.adjust.method = "BH", sort = TRUE): object 'pseq' not found
```


### Wilcoxon test (two-group comparisons)

If the data remarkably violates Gaussian assumptions, you need to use
non-parametric test. Wilcoxon is one option for two group
comparison. Here we compare males and females in the example data. 


```r
pval <- check_wilcoxon(pseq, "gender")
```

```
## Error in sample_data(x): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq' not found
```


