### ROC analysis

A basic example of ROC/AUC analysis with simulated random data.


```r
# Define two sample groups for demonstration purpose
g1 <- sample(colnames(oligo.data), 10)
g2 <- setdiff(colnames(oligo.data), g1)

# Compare the two groups with t-test
pvalues <- c()
for (tax in rownames(oligo.data)) {
  pvalues[[tax]] <- t.test(oligo.data[tax, g1], oligo.data[tax, g2])$p.value
}

# Order the taxa based on the p-values
ordered.results <- names(sort(pvalues))

# Assume there are some known true positives (here just randomly picked for demonstration)
true.positives <- sample(rownames(oligo.data), 10)

# Overall ROC analysis (this will give the cumulative TPR and FPR along the ordered list)
res <- roc(ordered.results, true.positives)

# Calculate ROC/AUC value
auc <- roc.auc(ordered.results, true.positives)
print(auc)
```

```
## [1] 0.6126459
```

```r
# Plot ROC curve
roc.plot(ordered.results, true.positives, line = TRUE)
```

![plot of chunk roc-example](figure/roc-example-1.png) 
