
### Correlation analysis

Cross-correlate input variables between two data sets (samples x features matrices; matched samples assumed; correlations calculated between features). For instance, between phylotypes and environmental variables:


```r
# Determine sample size
N <- ncol(l2.log10.simulated)

# Create some random sample information for the metadata
metadata.simulated <- data.frame(list(
              sampleID = colnames(l2.log10.simulated),
	      subjectID = paste("subjectID", rep(1:4, 5)),
	      group = sample(paste("group", rep(1:4, 5))),
	      time = paste("group", rep(1:4, 5)),
              age = runif(N, 0, 100),
              gender = sample(c("M", "F"), N, replace = TRUE),
              diet = sample(c("Apricots", "Beverages", "Carrots"), 
	      N, replace = TRUE)))
```

```
## Error in data.frame(sampleID = c("sample-1", "sample-2", "sample-3", "sample-4", : arguments imply differing number of rows: 44, 20
```

```r
# Calculate correlations between metadata and profiling data
cc <- cross.correlate(metadata.simulated, t(l2.log10.simulated), 
      			mode = "table", p.adj.method = "BY")
```

```
## Error in cor.test.default(xi, y[, j], method = method, use = "pairwise.complete.obs"): 'x' and 'y' must have the same length
```

```r
print(head(cc))
```

```
## Error in head(cc): error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'cc' not found
```



Matrix mode plus visualization for the similarity matrix:  


```r
cc <- cross.correlate(metadata.simulated, t(l2.log10.simulated), mode = "matrix")
```

```
## Error in cor.test.default(xi, y[, j], method = method, use = "pairwise.complete.obs"): 'x' and 'y' must have the same length
```

```r
tmp <- netresponse::plot_matrix(t(cc$cor), type = "twoway")
```

```
## Error in t(cc$cor): object 'cc' not found
```
