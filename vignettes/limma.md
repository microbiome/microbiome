### Limma analysis

Example of limma analysis. Identify most significantly different taxa between males and females. For further details, see [limma homepage](http://bioinf.wehi.edu.au/limma/) and [limma User's guide](http://www.lcg.unam.mx/~lcollado/R/resources/limma-usersguide.pdf). For discussion on why limma is preferred over t-test, see [this article](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Get example data
library(microbiome)
pseq <- download_microbiome("peerj32")$physeq
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

```r
otu <- log10(otu_table(pseq)@.Data)
meta <- sample_data(pseq)
groups <- meta$gender

# Compare the two groups with limma
library(limma)

# Prepare the design matrix which states the groups for each sample
# in the otu
design <- cbind(intercept = 1, Grp2vs1 = groups)
rownames(design) <- rownames(meta)
design <- design[colnames(otu), ]

# NOTE: results and p-values are given for all groupings in the design matrix
# Now focus on the second grouping ie. pairwise comparison
coef.index <- 2
     
# Fit the limma model
fit <- lmFit(otu, design)
fit <- eBayes(fit)

# Summarise or plot the results
topTable(fit, coef = coef.index)
```

```
##                                     logFC   AveExpr         t      P.Value
## Uncultured Clostridiales II    -0.4121655 1.3673693 -3.721676 0.0005453139
## Eubacterium siraeum et rel.    -0.3361407 1.6705368 -3.518974 0.0009997625
## Clostridium nexile et rel.      0.1821841 2.8405046  3.411676 0.0013694131
## Sutterella wadsworthia et rel. -0.3274127 1.4996997 -3.132421 0.0030375918
## Uncultured Clostridiales I     -0.6792917 1.4208331 -2.892738 0.0058532773
## Allistipes et rel.             -0.2514321 2.2670680 -2.874970 0.0061381919
## Aerococcus                      0.4209949 0.5788959  2.769463 0.0081131928
## Clostridium (sensu stricto)    -0.4319492 0.7262054 -2.696675 0.0098022788
## Eubacterium rectale et rel.     0.1365560 2.6448705  2.618545 0.0119712074
## Oxalobacter formigenes et rel. -0.2244230 2.1803707 -2.615543 0.0120627017
##                                 adj.P.Val          B
## Uncultured Clostridiales II    0.05934123 -0.2396518
## Eubacterium siraeum et rel.    0.05934123 -0.7698536
## Clostridium nexile et rel.     0.05934123 -1.0441290
## Sutterella wadsworthia et rel. 0.09872173 -1.7350274
## Uncultured Clostridiales I     0.13299416 -2.2988134
## Allistipes et rel.             0.13299416 -2.3394452
## Aerococcus                     0.15067358 -2.5772357
## Clostridium (sensu stricto)    0.15681512 -2.7377198
## Eubacterium rectale et rel.    0.15681512 -2.9066313
## Oxalobacter formigenes et rel. 0.15681512 -2.9130498
```

```r
# Q-Q plot
qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior)
abline(0,1)
```

![plot of chunk limma-example](figure/limma-example-1.png) 

```r
# Volcano plot
volcanoplot(fit, coef = coef.index, highlight = coef.index)
```

![plot of chunk limma-example](figure/limma-example-2.png) 

```r
# Adjusted p-values; show all significant ones
pvalues.limma <- p.adjust(fit$p.value[, coef.index], method = "fdr")
names(pvalues.limma) <- rownames(fit$p.value)
print(sort(pvalues.limma[pvalues.limma < 0.1]))
```

```
##     Clostridium nexile et rel.    Eubacterium siraeum et rel. 
##                     0.05934123                     0.05934123 
##    Uncultured Clostridiales II Sutterella wadsworthia et rel. 
##                     0.05934123                     0.09872173
```


### Comparison between limma and t-test

Order the taxa with t-test for comparison and validation purposes. The
differences are small in this simulated example, but [can be
considerable in real
data](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Compare the two groups with t-test
pvalues.ttest <- c()
male.samples <- filter(meta, gender == "male")$sample
```

```
## Error in filter(meta, gender == "male"): object 'gender' not found
```

```r
female.samples <- filter(meta, gender == "female")$sample
```

```
## Error in filter(meta, gender == "female"): object 'gender' not found
```

```r
for (tax in rownames(otu)) {
  pvalues.ttest[[tax]] <- t.test(otu[tax, male.samples], otu[tax, female.samples])$p.value
}
```

```
## Error in t.test(otu[tax, male.samples], otu[tax, female.samples]): object 'male.samples' not found
```

```r
# Multiple testing correction
pvalues.ttest <- p.adjust(pvalues.ttest, method = "fdr")

# Compare p-values between limma and t-test
taxa <- rownames(otu)
plot(pvalues.ttest[taxa], pvalues.limma[taxa])
```

```
## Error in plot.window(...): need finite 'xlim' values
```

![plot of chunk limma-compairson](figure/limma-compairson-1.png) 

```r
abline(0,1,lty = 2)
```

### TODO

Check relations to expressionSets. Could we use tools directly from that context by a suitable conversion.
