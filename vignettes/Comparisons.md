

## Group-wise comparisons


Read example data


```r
# Example data
library(microbiome)
data(peerj32)
x <- peerj32$microbes
m <- peerj32$meta
```

### Comparing of two or more groups with a parametric test (linear model; ANOVA)


```r
# Define sample groups (gender + case/ctrl)
group <- apply(cbind(m$gender, m$group), 1, function (x) {paste(x, collapse = "-")})
names(group) <- m$sample

# Calculate p-values for multi-group comparison
# NOTE: the HITChip data (x) must be at log10 scale in parametric tests!
x <- log10(t(x))
pvals <- check.anova(x, group, p.adjust.method = "BH", sort = TRUE)

# Writing the results to a table
# write.table(pvals, file = "myfile.tab", quote = F)
```


### Wilcoxon test (only for two-group comparison)

If the data remarkably violates Gaussian assumptions, you need to use
non-parametric test, Wilcoxon is one option for two group
comparison. Here we compare males and females in simulated toy
data. The sampleIDs in G1 and G2 must match with input matrix
columns. Compare two sample groups and return adjusted p-values for
each phylotype. This test uses unpaired Wilcoxon test with BH
correction.


```r
# Calculate Wilcoxon test with BH p-value correction for gender
pval <- apply(x, 1, function (xi) {wilcox.test(xi ~ m$gender)$p.value})

# Multiple testing correction
pval <- p.adjust(pval, "fdr")

# Sort p-values
pval <- sort(pval)

# Investigate the top findings
head(pval)
```

```
##           Clostridium nexile et rel.         Megasphaera elsdenii et rel. 
##                           0.04177601                           0.04177601 
##          Uncultured Clostridiales II          Eubacterium siraeum et rel. 
##                           0.04177601                           0.08356243 
##        Dorea formicigenerans et rel. Outgrouping clostridium cluster XIVa 
##                           0.08797056                           0.09743840
```


