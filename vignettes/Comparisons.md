
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
group <- meta$group
names(group) <- meta$sample

# Calculate 1-way ANOVA p-values for the multi-group comparison
pvals <- check_anova(pseq, group, p.adjust.method = "BH", sort = TRUE)
```


### Wilcoxon test (two-group comparisons)

If the data remarkably violates Gaussian assumptions, you need to use
non-parametric test. Wilcoxon is one option for two group
comparison. Here we compare males and females in the example data. 


```r
pval <- check_wilcoxon(pseq, "gender")
```


