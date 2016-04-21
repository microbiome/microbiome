### Limma analysis

Example of limma analysis. Identify most significantly different taxa between males and females. For further details, see [limma homepage](http://bioinf.wehi.edu.au/limma/) and [limma User's guide](http://www.lcg.unam.mx/~lcollado/R/resources/limma-usersguide.pdf). For discussion on why limma is preferred over t-test, see [this article](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Get example data
library(microbiome)
data("peerj32")
pseq <- peerj32$phyloseq
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

# Summarise 
kable(topTable(fit, coef = coef.index), digits = 2)
```



|                               | logFC| AveExpr|     t| P.Value| adj.P.Val|     B|
|:------------------------------|-----:|-------:|-----:|-------:|---------:|-----:|
|Uncultured Clostridiales II    | -0.41|    1.37| -3.72|    0.00|      0.06| -0.24|
|Eubacterium siraeum et rel.    | -0.34|    1.67| -3.52|    0.00|      0.06| -0.77|
|Clostridium nexile et rel.     |  0.18|    2.84|  3.41|    0.00|      0.06| -1.04|
|Sutterella wadsworthia et rel. | -0.33|    1.50| -3.13|    0.00|      0.10| -1.74|
|Uncultured Clostridiales I     | -0.68|    1.42| -2.89|    0.01|      0.13| -2.30|
|Allistipes et rel.             | -0.25|    2.27| -2.87|    0.01|      0.13| -2.34|
|Aerococcus                     |  0.42|    0.58|  2.77|    0.01|      0.15| -2.58|
|Clostridium (sensu stricto)    | -0.43|    0.73| -2.70|    0.01|      0.16| -2.74|
|Eubacterium rectale et rel.    |  0.14|    2.64|  2.62|    0.01|      0.16| -2.91|
|Oxalobacter formigenes et rel. | -0.22|    2.18| -2.62|    0.01|      0.16| -2.91|


Adjusted p-values; show all significant ones:


```r
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


### Q-Q plot



```r
qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior)
abline(0,1)
```

![plot of chunk limma-qq](figure/limma-qq-1.png)

### Volcano plot


```r
volcanoplot(fit, coef = coef.index, highlight = coef.index)
```

![plot of chunk limma-volcano](figure/limma-volcano-1.png)



### Comparison between limma and t-test

Order the taxa with t-test for comparison and validation purposes. The
differences are small in this simulated example, but [can be
considerable in real
data](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Compare the two groups with t-test
library(dplyr)
pvalues.ttest <- c()
male.samples <- dplyr::filter(meta, gender == "male")$sample
female.samples <- dplyr::filter(meta, gender == "female")$sample
for (tax in rownames(otu)) {
  pvalues.ttest[[tax]] <- t.test(otu[tax, male.samples], otu[tax, female.samples])$p.value
}
# Multiple testing correction
pvalues.ttest <- p.adjust(pvalues.ttest, method = "fdr")

# Compare p-values between limma and t-test
taxa <- rownames(otu)
plot(pvalues.ttest[taxa], pvalues.limma[taxa])
abline(0,1,lty = 2)
```

![plot of chunk limma-compairson](figure/limma-compairson-1.png)

### Continuous variables

Quantify continuous associations with lm_phyloseq function. This uses
the limma model to generate a table of P-values and effect sizes. Note
that no confounding variables taken into account in this wrapper. See
the [limma homepage](http://bioinf.wehi.edu.au/limma/) for more
detailed analyses.


```r
data("atlas1006")
tab <- lm_phyloseq(atlas1006, "age")
```

```
## Error in qr.default(x): NA/NaN/Inf in foreign function call (arg 1)
```

```r
kable(head(tab))
```

```
## Error in head(tab): error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'tab' not found
```


