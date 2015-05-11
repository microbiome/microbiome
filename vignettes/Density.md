### Abundance histograms

Different sample sets have different population distributions in
microbial abundance. It is also important to consider whether to use
absolute or logarithmic abundances!

Start by loading the [HITChip Atlas data set](Data.Rmd)



```r
# Load microbiome package
library(microbiome)  

# Load tools
library(dplyr)

# 1. List all samples (all time points and DNA extraction methods)
all.samples <- meta$SampleID

# 2. List samples at time point 0 that have specific DNA extraction method 
rbb.samples <- filter(meta, Time == "0" & DNA_extraction_method == "r")$SampleID
```

```
## Error in filter_impl(.data, dots): object 'Time' not found
```

```r
# Visualize
#tax <- "Prevotella.melaninogenica.et.rel."
tax <- "Bifidobacterium"
d <- data[all.samples, tax]
par(mfrow = c(1, 2))
plot(density(d), main = paste(tax, "(All samples)"), xlab = "Abundance (Absolute HITChip signal)")
```

```
## Error in plot(density(d), main = paste(tax, "(All samples)"), xlab = "Abundance (Absolute HITChip signal)"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error in density.default(d) : 
##   need at least 2 points to select a bandwidth automatically
```

```r
plot(density(log10(d)), main = paste(tax, "(All samples)"), xlab = "Abundance (Log10 HITChip signal)")
```

```
## Error in plot(density(log10(d)), main = paste(tax, "(All samples)"), xlab = "Abundance (Log10 HITChip signal)"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error in density.default(log10(d)) : 
##   need at least 2 points to select a bandwidth automatically
```

```r
d <- data[rbb.samples, tax]
```

```
## Error in eval(expr, envir, enclos): object 'rbb.samples' not found
```

```r
par(mfrow = c(1, 2))
plot(density(d), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Absolute HITChip signal)")
```

```
## Error in plot(density(d), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Absolute HITChip signal)"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error in density.default(d) : 
##   need at least 2 points to select a bandwidth automatically
```

```r
plot(density(log10(d)), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Log10 HITChip signal)")
```

```
## Error in plot(density(log10(d)), main = paste(tax, "(RBB samples)"), xlab = "Abundance (Log10 HITChip signal)"): error in evaluating the argument 'x' in selecting a method for function 'plot': Error in density.default(log10(d)) : 
##   need at least 2 points to select a bandwidth automatically
```
