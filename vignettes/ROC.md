### ROC analysis

A basic example of ROC/AUC analysis with simulated random data.


```r
# Define two sample groups for demonstration purpose
g1 <- sample(colnames(oligo.data), 10)
```

```
## Error in colnames(oligo.data): error in evaluating the argument 'x' in selecting a method for function 'colnames': Error: object 'oligo.data' not found
```

```r
g2 <- setdiff(colnames(oligo.data), g1)
```

```
## Error in setdiff(colnames(oligo.data), g1): error in evaluating the argument 'x' in selecting a method for function 'setdiff': Error in colnames(oligo.data) : 
##   error in evaluating the argument 'x' in selecting a method for function 'colnames': Error: object 'oligo.data' not found
```

```r
# Compare the two groups with t-test
pvalues <- c()
for (tax in rownames(oligo.data)) {
  pvalues[[tax]] <- t.test(oligo.data[tax, g1], oligo.data[tax, g2])$p.value
}
```

```
## Error in rownames(oligo.data): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'oligo.data' not found
```

```r
# Order the taxa based on the p-values
ordered.results <- names(sort(pvalues))

# Assume there are some known true positives (here just randomly picked for demonstration)
true.positives <- sample(rownames(oligo.data), 10)
```

```
## Error in rownames(oligo.data): error in evaluating the argument 'x' in selecting a method for function 'rownames': Error: object 'oligo.data' not found
```

```r
# Overall ROC analysis (this will give the cumulative TPR and FPR along the ordered list)
res <- roc(ordered.results, true.positives)
```

```
## Error in roc(ordered.results, true.positives): object 'true.positives' not found
```

```r
# Calculate ROC/AUC value
auc <- roc.auc(ordered.results, true.positives)
```

```
## Error in roc(ordered.results, true.positives): object 'true.positives' not found
```

```r
print(auc)
```

```
## Error in print(auc): object 'auc' not found
```

```r
# Plot ROC curve
roc.plot(ordered.results, true.positives, line = TRUE)
```

```
## Error in roc(ordered.results, true.positives): object 'true.positives' not found
```
