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
##                                  logFC  AveExpr          t      P.Value
## Methylobacterium             -1.958649 2.635062 -50.137611 3.246868e-22
## Aerococcus                   -1.964847 2.670740 -47.823533 8.160712e-22
## Bulleidia moorei et rel.     -2.014792 3.350158 -45.135816 2.519169e-21
## Bacillus                     -1.937091 3.176410 -41.287269 1.427204e-20
## Anaerobiospirillum           -1.973006 2.657116 -40.663333 1.918863e-20
## Haemophilus                  -1.972022 2.987658 -30.163588 6.222907e-18
## Actinomycetaceae             -1.937329 3.006550 -26.722889 6.373711e-17
## Prevotella tannerae et rel.  -2.003901 3.707903 -10.208120 2.669446e-09
## Bacteroides plebeius et rel. -1.942524 3.907811 -10.088794 3.246699e-09
## Uncultured Clostridiales I   -1.997527 4.515231  -7.413681 4.122638e-07
##                                 adj.P.Val         B
## Methylobacterium             4.220928e-20 41.162031
## Aerococcus                   5.304463e-20 40.238565
## Bulleidia moorei et rel.     1.091640e-19 39.104023
## Bacillus                     4.638411e-19 37.348552
## Anaerobiospirillum           4.989044e-19 37.047867
## Haemophilus                  1.348296e-16 31.126439
## Actinomycetaceae             1.183689e-15 28.724099
## Prevotella tannerae et rel.  4.337850e-08 10.474770
## Bacteroides plebeius et rel. 4.689676e-08 10.271314
## Uncultured Clostridiales I   5.359430e-06  5.250828
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
##  [1] 4.220928e-20 5.304463e-20 1.091640e-19 4.638411e-19 4.989044e-19
##  [6] 1.348296e-16 1.183689e-15 4.337850e-08 4.689676e-08 5.359430e-06
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
