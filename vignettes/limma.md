### Limma analysis

Example of limma analysis with simulated random data. For further details, see [limma homepage](http://bioinf.wehi.edu.au/limma/) and [limma User's guide](http://www.lcg.unam.mx/~lcollado/R/resources/limma-usersguide.pdf). For discussion on why limma is preferred over t-test, see [this article](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).



```r
# Get example data
library(microbiome, quietly = TRUE)

# Define here your own HITChip data folder
#data.directory <- system.file("extdata", package = "microbiome")
data.directory <- "~/R/x86_64-pc-linux-gnu-library/3.0/microbiome/extdata"

# Read HITChip data
hitchip.matrix <- read.profiling(level = "L2", 
                       data.dir = data.directory, log10 = TRUE)
```

```
## Reading ~/R/x86_64-pc-linux-gnu-library/3.0/microbiome/extdata/L2-frpa.tab
## Logarithmizing the data
```

```r
# Define two random groups for demonstration purpose
g1 <- sample(colnames(hitchip.matrix), 10)
g2 <- setdiff(colnames(hitchip.matrix), g1)
# Modify hitchip matrix so that there are a few significant differences
altered.taxa <- sample(rownames(hitchip.matrix), 10)
hitchip.matrix[altered.taxa, g1] <- hitchip.matrix[altered.taxa, g1] + 2

# Compare the two groups with limma
library(limma)
```

```
## 
## Attaching package: 'limma'
## 
## The following object is masked from 'package:BiocGenerics':
## 
##     plotMA
```

```r
# Prepare the design matrix which states the groups for each sample
# in the hitchip.matrix
design <- cbind(intercept=1, Grp2vs1=c(rep(0, length(g1)), rep(1, length(g2))))
rownames(design) <- c(g1, g2)
design <- design[colnames(hitchip.matrix), ]

# NOTE: results and p-values are given for all groupings in the design matrix
# Now focus on the second grouping ie. pairwise comparison
coef.index <- 2
     
# Fit the limma model
fit <- lmFit(hitchip.matrix, design)
fit <- eBayes(fit)

# Summarise or plot the results
topTable(fit, coef = coef.index)
```

```
##                                      logFC  AveExpr          t
## Eubacterium limosum et rel.      -1.996643 3.178103 -44.244703
## Wissella et rel.                 -2.087264 2.758124 -29.878412
## Lactococcus                      -2.048501 3.177703 -29.824249
## Lactobacillus plantarum et rel.  -2.125528 3.963640 -27.460242
## Bryantella formatexigens et rel. -2.243230 4.517741 -13.750185
## Uncultured Clostridiales II      -2.275919 4.440935 -11.383575
## Uncultured Mollicutes            -1.935851 3.977371 -10.831144
## Clostridium ramosum et rel.      -1.732479 3.736636  -9.550878
## Ruminococcus obeum et rel.       -2.044862 5.535052  -9.358625
## Streptococcus mitis et rel.      -1.935691 4.358083  -8.697516
##                                       P.Value    adj.P.Val         B
## Eubacterium limosum et rel.      3.378865e-21 4.392524e-19 38.806959
## Wissella et rel.                 6.931890e-18 3.110880e-16 31.018845
## Lactococcus                      7.178954e-18 3.110880e-16 30.982758
## Lactobacillus plantarum et rel.  3.524020e-17 1.145306e-15 29.341105
## Bryantella formatexigens et rel. 1.471547e-11 3.826021e-10 15.895750
## Uncultured Clostridiales II      4.084953e-10 8.850731e-09 12.436140
## Uncultured Mollicutes            9.570763e-10 1.777427e-08 11.550594
## Clostridium ramosum et rel.      7.814506e-09 1.269857e-07  9.368664
## Ruminococcus obeum et rel.       1.088754e-08 1.572645e-07  9.024415
## Streptococcus mitis et rel.      3.525558e-08 4.583225e-07  7.805693
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
names(pvalues.limma) <- unname(unlist(fit$genes))
print(sort(pvalues.limma[pvalues.limma < 0.05]))
```

```
##  [1] 4.392524e-19 3.110880e-16 3.110880e-16 1.145306e-15 3.826021e-10
##  [6] 8.850731e-09 1.777427e-08 1.269857e-07 1.572645e-07 4.583225e-07
```

### Comparison between limma and t-test

Order the taxa with t-test for comparison and validation purposes. The
differences are small in this simulated example, but [can be
considerable in real
data](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).


```r
# Compare the two groups with t-test
pvalues.ttest <- c()
for (tax in rownames(hitchip.matrix)) {
  pvalues.ttest[[tax]] <- t.test(hitchip.matrix[tax, g1], hitchip.matrix[tax, g2])$p.value
}
# Multiple testing correction
pvalues.ttest <- p.adjust(pvalues.ttest, method = "fdr")


# Order the taxa based on the p-values
taxa <- rownames(hitchip.matrix)
plot(pvalues.ttest[taxa], pvalues.limma[taxa])
```

```
## Error in plot.window(...): need finite 'ylim' values
```

![plot of chunk limma-compairson](figure/limma-compairson-1.png) 

```r
abline(0,1,lty = 2)
```
