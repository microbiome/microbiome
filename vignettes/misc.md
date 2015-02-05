### Writing diversity table into file


```r
output.dir <- "./"
write.table(div.table, file = "DiversityTable.tab", sep = "\t")
```

## Save clustering image to a file

Save in PDF:


```r
pdf("myplot.pdf", width = 7, height = 7 * length(hc$order)/20)
plot(hc, hang=-1, main = "Hierarchical clustering")
dev.off()
```

Save in TIFF:


```r
tiff("myplot.tif", width = 480, height = 480 * length(hc$order)/20)
plot(hc, hang=-1, main = "Hierarchical clustering")
dev.off()
```

To save in Microsoft EMF format, try the following. If you find a
way to tune figure width for emf files kindly let the admins know.


```r
plot(hc, hang=-1, main = "Hierarchical clustering")
savePlot("myplot.emf", type = "emf")
dev.off()
```

## Cluster significance testing

Cluster significance can be assessed with bootstrap resampling of the
original data. Here, perform R = 1000 bootstrap iterations and return
significant classes with at least 2 members and correlation higher
than specified. Only significant clusters with the given p-value
threshold are returned; p-values are corrected for multiple testing
with Benjamini-Hochberg method.


```r
significant.clusters <- hclust.significance(t(mydata), R = 1e3, min.size = 2, corr.th = 0.5, metric = "pearson", pvalue.threshold = 0.05) 
```

```
## Error in eval(expr, envir, enclos): could not find function "hclust.significance"
```

