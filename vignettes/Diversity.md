
To run these examples, you need to download the [HITChip Atlas data set](Data.md)

### Richness 


```r
# Use the HITChip Atlas data
mydata <- atlas$microbes
```

```
## Error in eval(expr, envir, enclos): object 'atlas' not found
```

```r
# Determine detection threshold at the 0.15 quantile
# Bacteria that exceed this threshold are considered present
# otherwise absent
det.th <- quantile(mydata, 0.15)

# Visualize the detection threshold (at log10 scale for clarity)
plot(density(log10(mydata)), main = "Detection threshold", xlab = "Abundance (Log10)", ylab = "Frequency")
abline(v = log10(det.th))
```

![plot of chunk rich-example](figure/rich-example-1.png) 

```r
# Calculate richness. 
# This indicates how many oligos are present in each sample
# (exceed the detection threshold)
ri <- rowSums(mydata > det.th)
hist(ri, main = "Richness")
```

![plot of chunk rich-example](figure/rich-example-2.png) 


### Diversity 


```r
# Diversity using the vegan package
# NOTE: data needs to be in absolute scale, not logarithmic
di <- vegan::diversity(mydata, index = "shannon")
hist(di, main = "Diversity")
```

![plot of chunk div-example](figure/div-example-1.png) 

### Evenness


```r
# Pielou's evenness is S/ln(R) w.r.t. given detection threshold
# NOTE: here we use detection threshold for diversity as well because
# the exact same data has to be used for diversity and richness calculation,
# and for richness calculation the detection threshold needs to be set anyway
# Diversity can be as such calculated also without threshold (see above)
# but that gives somewhat different result.
mydata2 <- mydata - det.th # NOTE: absolute (not log) scale data
S <- vegan::diversity(mydata2, index = "shannon")
R <- rowSums(mydata2 > 0)
ev <- S/log(R)
```


### Compare with known background factors

Diversity vs. obesity


```r
par(mar = c(6, 4, 3, 1))
bmi <- atlas$meta$BMI_group
```

```
## Error in eval(expr, envir, enclos): object 'atlas' not found
```

```r
di <- vegan::diversity(atlas$microbes)
```

```
## Error in as.matrix(x): object 'atlas' not found
```

```r
boxplot(di ~ bmi, las = 2, main = "Microbiota diversity vs. obesity")
```

```
## Error in eval(expr, envir, enclos): object 'bmi' not found
```


Diversity vs. age with smoothed confidence intervals:


```r
library(microbiome)
library(sorvi)
library(dplyr)

# Pick the subset of RBB-preprocessed samples from time point 0
rbb.samples <- filter(atlas$meta, Time == 0 & DNA_extraction_method == "r")$SampleID
```

```
## Error in filter_(.data, .dots = lazyeval::lazy_dots(...)): object 'atlas' not found
```

```r
# Collect variables into a data frame
df <- data.frame(Age = atlas$meta[rbb.samples, "Age"], Diversity = di[rbb.samples])
```

```
## Error in data.frame(Age = atlas$meta[rbb.samples, "Age"], Diversity = di[rbb.samples]): object 'atlas' not found
```

```r
# Visualize
p <- sorvi::regression_plot(Diversity~Age, df, shade = TRUE, mweight = TRUE, verbose = FALSE)
```

```
## Error in eval(expr, envir, enclos): incorrect size (1), expecting : 44
```

```r
print(p)
```

![plot of chunk visu-example3](figure/visu-example3-1.png) 

