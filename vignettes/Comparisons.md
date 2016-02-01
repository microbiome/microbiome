
## Group-wise comparisons

Read example data from a [diet swap study](http://dx.doi.org/10.1038/ncomms7342):


```r
library(microbiome)
data("dietswap")
pseq <- dietswap
```

### Comparing of two or more groups with a parametric test (linear model; ANOVA)


```r
# 1-way ANOVA p-values for the multi-group comparison across time groups
anova.results <- check_anova(pseq, "group", p.adjust.method = "BH")
kable(head(anova.results))
```



|   p.anova|   p.ED-DI|   p.HE-DI|   p.HE-ED|      ave.DI|      ave.ED|    ave.HE|
|---------:|---------:|---------:|---------:|-----------:|-----------:|---------:|
| 0.0000095| 0.9575836| 0.0000040| 0.0000008| 10942.94000| 10251.71000| 23198.120|
| 0.0012025| 0.7742956| 0.0000631| 0.0007913|  3552.40300|  4366.57300|  8740.280|
| 0.0012025| 0.8718784| 0.0000869| 0.0005321|  1344.44400|  1495.02700|  2632.413|
| 0.0016075| 0.8299227| 0.0011873| 0.0001083|  4252.19400|  3543.70700|  8621.707|
| 0.0076210| 0.7415585| 0.0004481| 0.0052958|   944.22220|  1109.60000|  1808.480|
| 0.0088024| 0.0119282| 0.0003993| 0.5721161|    61.02778|    66.65333|    68.600|


### Wilcoxon test (two-group comparisons)

If the data remarkably violates Gaussian assumptions use
non-parametric test. Wilcoxon is one option for two group
comparison. Here we compare males and females in the example data.


```r
pval <- check_wilcoxon(pseq, "sex")
```



