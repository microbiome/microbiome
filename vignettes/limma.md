## Linear models with limma


### Discrete variables: sex

Identify most significantly different taxa between males and females.

For further details, see [limma
homepage](http://bioinf.wehi.edu.au/limma/) and [limma User's
guide](http://www.lcg.unam.mx/~lcollado/R/resources/limma-usersguide.pdf). For
discussion on why limma is preferred over t-test, see [this
article](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Get example data
library(microbiome)
data("peerj32")
pseq <- peerj32$phyloseq
otu <- transform_phyloseq(get_taxa(pseq), "log10")
```

```
## Error in .nextMethod(.Object = .Object, ... = ...): argument "taxa_are_rows" is missing, with no default
```

```r
meta <- sample_data(pseq)
groups <- meta$gender

# Compare the two groups with limma
library(limma)

# Prepare the design matrix which states the groups for each sample
# in the otu
design <- cbind(intercept = 1, Grp2vs1 = groups)
rownames(design) <- rownames(meta)
design <- design[colnames(otu), ]
```

```
## Error in design[colnames(otu), ]: subscript out of bounds
```

```r
# NOTE: results and p-values are given for all groupings in the design matrix
# Now focus on the second grouping ie. pairwise comparison
coef.index <- 2
     
# Fit the limma model
fit <- lmFit(otu, design)
```

```
## Error in lmFit(otu, design): row dimension of design doesn't match column dimension of data object
```

```r
fit <- eBayes(fit)
```

```
## Error in ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, : object 'fit' not found
```

```r
# Summarise 
kable(topTable(fit, coef = coef.index, p.value=0.05), digits = 2)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```


### Q-Q plot



```r
qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior)
```

```
## Error in qqt(fit$t[, coef.index], df = fit$df.residual + fit$df.prior): object 'fit' not found
```

```r
abline(0,1)
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

### Volcano plot


```r
volcanoplot(fit, coef = coef.index, highlight = coef.index)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```



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
```

```
## Error in otu[tax, male.samples]: subscript out of bounds
```

```r
# Multiple testing correction
pvalues.ttest <- p.adjust(pvalues.ttest, method = "fdr")

# Compare p-values between limma and t-test
taxa <- rownames(otu)
plot(pvalues.ttest[taxa], pvalues.limma[taxa])
```

```
## Error in plot(pvalues.ttest[taxa], pvalues.limma[taxa]): error in evaluating the argument 'y' in selecting a method for function 'plot': Error: object 'pvalues.limma' not found
```

```r
abline(0,1,lty = 2)
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

### Continuous variables

Quantify continuous associations with lm_phyloseq function. This uses
the limma model to generate a table of P-values and effect sizes. Note
that no confounding variables taken into account in this wrapper. See
the [limma homepage](http://bioinf.wehi.edu.au/limma/) for more
detailed analyses.


```r
data("atlas1006")
source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
tab <- lm_phyloseq(atlas1006, "age")
kable(head(tab), digits = 3)
```



|                                   |  logFC| AveExpr|       t| P.Value| adj.P.Val|      B|
|:----------------------------------|------:|-------:|-------:|-------:|---------:|------:|
|Bifidobacterium                    | -0.015|   3.702| -12.507|       0|         0| 63.502|
|Clostridium difficile et rel.      | -0.009|   3.229|  -9.890|       0|         0| 37.157|
|Oscillospira guillermondii et rel. |  0.012|   4.535|   9.828|       0|         0| 36.594|
|Bacteroides splachnicus et rel.    |  0.006|   3.219|   9.552|       0|         0| 34.132|
|Collinsella                        | -0.009|   2.828|  -9.107|       0|         0| 30.273|
|Tannerella et rel.                 |  0.007|   3.162|   8.977|       0|         0| 29.175|


