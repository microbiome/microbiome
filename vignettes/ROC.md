### ROC analysis

A basic example of ROC/AUC analysis.


### Load example data


```r
library(microbiome)
#pseq <- download_microbiome("dietswap")
data("dietswap")
pseq <- dietswap

# Pick two groups of samples
# African, DI group, time points 1 and 2
# See the original publication for details: 
# references provided in help(dietswap)
pseq1 <- subset_samples(pseq, nationality == "AFR" & 
		     timepoint.within.group == 1 & 
		     group == "DI")
pseq2 <- subset_samples(pseq, nationality == "AFR" & 
		     timepoint.within.group == 2 & 
		     group == "DI")

# Pick OTU matrix
otu <- otu_table(pseq)@.Data

# Pick sample metadata
meta <- sample_data(pseq)

# Define two sample groups for demonstration purpose
g1 <- sample_names(pseq1)
g2 <- sample_names(pseq2)

# Compare the two groups with Wilcoxon test
pvalues <- c()
for (tax in rownames(otu)) {
  pvalues[[tax]] <- wilcox.test(otu[tax, g1], otu[tax, g2])$p.value
}

# Assume there are some known true positives 
# Here for instance Bacteroidetes
bacteroidetes <- levelmap("Bacteroidetes", from = "L1", to = "L2", GetPhylogeny("HITChip", "filtered"))$Bacteroidetes
```


### Overall ROC analysis 

Based on the [xrobin/pROC](https://github.com/xrobin/pROC) package
(see that page for more examples with confidence limits etc):


```r
library(pROC)
```

```
## Type 'citation("pROC")' for a citation.
```

```
## 
## Attaching package: 'pROC'
```

```
## The following objects are masked from 'package:stats':
## 
##     cov, smooth, var
```

```r
res <- roc(names(pvalues) %in% bacteroidetes, pvalues)
```


### ROC/AUC value


```r
res$auc
```

```
## Area under the curve: 0.6373
```


### Plot ROC curve


```r
plot(res)
```

![plot of chunk roc-example4](figure/roc-example4-1.png)

```
## 
## Call:
## roc.default(response = names(pvalues) %in% bacteroidetes, predictor = pvalues)
## 
## Data: pvalues in 107 controls (names(pvalues) %in% bacteroidetes FALSE) < 16 cases (names(pvalues) %in% bacteroidetes TRUE).
## Area under the curve: 0.6373
```
