### Interactive motioncharts with googleVis

This provides interactive charts where you can change the variables at
X/Y axis, point size and color as you like. With time series data
moving animations are also possible.


Prepare the data:


```r
library(microbiome)  

data(peerj32)

# Species-level data in phyloseq format
pseq <- peerj32$phyloseq

# Retrieve L2 (genus-like) data in phyloseq format
pseq.L2 <- aggregate_taxa(pseq, level = "L2")
```

```
## Error in tax_glom(pseq, level): Bad taxrank argument. Must be among the values of rank_names(physeq)
```

```r
# Convert L2 to matrix format
genus.matrix.log10.simulated <- log10(otu_table(pseq.L2)@.Data)
```

```
## Error in otu_table(pseq.L2): error in evaluating the argument 'object' in selecting a method for function 'otu_table': Error: object 'pseq.L2' not found
```

```r
# Combine phylotype profiling data and sample metadata
df <- cbind(metadata.simulated, t(genus.matrix.log10.simulated))  
```

```
## Error in eval(expr, envir, enclos): object 'metadata.simulated' not found
```

Plot a Motion Chart using googleVis - package. Note: this requires
flash and internet connection. 

Form a motion chart from example data NOTE: the data set must be given
as data.frame which can contain NUMERIC and CHARACTER fields (NO
FACTORS, NOR LOGICAL variables!). The FIRST FOUR FIELDS must be
provided in the following order: idvar, timevar, two numeric fields,
then any number of numeric and character fields.

The plot shows only the first time point.  Replace the time field with
a constant to plot all in one figure using df$time <- rep(1, nrow(df))


```r
library(googleVis)  

# See help(gvisMotionChart) for further details
mchart <- gvisMotionChart(df, idvar = "sample", timevar = "time")  

# Plot immediately (opens in browser, requires flash)
plot(mchart)  
```

Save as html (needs javascript to open!). NOTE: html file viewing does not work locally - store the html file on server and view through internet:


```r
print(mchart, file="~/MotionChart.html")
```
