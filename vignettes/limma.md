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
```

```
## Error in eval(expr, envir, enclos): could not find function "otu_table"
```

```r
meta <- sample_data(pseq)
```

```
## Error in eval(expr, envir, enclos): could not find function "sample_data"
```

```r
groups <- meta$gender
```

```
## Error in eval(expr, envir, enclos): object 'meta' not found
```

```r
# Compare the two groups with limma
library(limma)
```

```
## 
## Attaching package: 'limma'
```

```
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
# Prepare the design matrix which states the groups for each sample
# in the otu
design <- cbind(intercept = 1, Grp2vs1 = groups)
```

```
## Error in eval(expr, envir, enclos): object 'groups' not found
```

```r
rownames(design) <- rownames(meta)
```

```
## Error in rownames(meta): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'meta' not found
```

```r
design <- design[colnames(otu), ]
```

```
## Error in design[colnames(otu), ]: error in evaluating the argument 'i' in selecting a method for function '[': Error in colnames(otu) : 
##   error in evaluating the argument 'x' in selecting a method for function 'colnames': Error: object 'otu' not found
```

```r
# NOTE: results and p-values are given for all groupings in the design matrix
# Now focus on the second grouping ie. pairwise comparison
coef.index <- 2
     
# Fit the limma model
fit <- lmFit(otu, design)
```

```
## Error in is(object, "list"): object 'otu' not found
```

```r
fit <- eBayes(fit)
```

```
## Error in ebayes(fit = fit, proportion = proportion, stdev.coef.lim = stdev.coef.lim, : object 'fit' not found
```

```r
# Summarise or plot the results
topTable(fit, coef = coef.index)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```

```r
# Q-Q plot
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

```r
# Volcano plot
volcanoplot(fit, coef = coef.index, highlight = coef.index)
```

```
## Error in is(fit, "MArrayLM"): object 'fit' not found
```

```r
# Adjusted p-values; show all significant ones
pvalues.limma <- p.adjust(fit$p.value[, coef.index], method = "fdr")
```

```
## Error in p.adjust(fit$p.value[, coef.index], method = "fdr"): object 'fit' not found
```

```r
names(pvalues.limma) <- rownames(fit$p.value)
```

```
## Error in rownames(fit$p.value): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'fit' not found
```

```r
print(sort(pvalues.limma[pvalues.limma < 0.1]))
```

```
## Error in sort(pvalues.limma[pvalues.limma < 0.1]): error in evaluating the argument 'x' in selecting a method for function 'sort': Error: object 'pvalues.limma' not found
```


### Comparison between limma and t-test

Order the taxa with t-test for comparison and validation purposes. The
differences are small in this simulated example, but [can be
considerable in real
data](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Compare the two groups with t-test
library(dplyr)
```

```
## 
## Attaching package: 'dplyr'
```

```
## The following object is masked from 'package:Biobase':
## 
##     combine
```

```
## The following objects are masked from 'package:BiocGenerics':
## 
##     combine, intersect, setdiff, union
```

```
## The following objects are masked from 'package:stats':
## 
##     filter, lag
```

```
## The following objects are masked from 'package:base':
## 
##     intersect, setdiff, setequal, union
```

```r
pvalues.ttest <- c()
male.samples <- dplyr::filter(meta, gender == "male")$sample
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'meta' not found
```

```r
female.samples <- dplyr::filter(meta, gender == "female")$sample
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'meta' not found
```

```r
for (tax in rownames(otu)) {
  pvalues.ttest[[tax]] <- t.test(otu[tax, male.samples], otu[tax, female.samples])$p.value
}
```

```
## Error in rownames(otu): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'otu' not found
```

```r
# Multiple testing correction
pvalues.ttest <- p.adjust(pvalues.ttest, method = "fdr")

# Compare p-values between limma and t-test
taxa <- rownames(otu)
```

```
## Error in rownames(otu): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'otu' not found
```

```r
plot(pvalues.ttest[taxa], pvalues.limma[taxa])
```

```
## Error in plot(pvalues.ttest[taxa], pvalues.limma[taxa]): object 'taxa' not found
```

```r
abline(0,1,lty = 2)
```

```
## Error in int_abline(a = a, b = b, h = h, v = v, untf = untf, ...): plot.new has not been called yet
```

### TODO

Check relations to expressionSets. Could we use tools directly from that context by a suitable conversion.
