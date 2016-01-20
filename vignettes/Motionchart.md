### Interactive motioncharts with googleVis

This provides interactive charts where you can change the variables at
X/Y axis, point size and color as you like. With time series data
moving animations are also possible.


Prepare the data:


```r
library(microbiome)  
```

```
## Warning: replacing previous import by 'ggplot2::Position' when loading
## 'phyloseq'
```

```
## Warning: replacing previous import by 'scales::alpha' when loading
## 'phyloseq'
```

```r
# Folder containing simulated example data
# (change your own path here)
data.directory <- system.file("extdata", package = "microbiome")

# Read HITChip profiling results
res <- read_hitchip(data.directory, method = "frpa")

# Pick the metadata (the meta.tab file in HITChip data folder)
metadata.simulated <- sample_data(res$pseq)
```

```
## Error in eval(expr, envir, enclos): could not find function "sample_data"
```

```r
# Species-level data in phyloseq format
pseq <- res$pseq 

# Retrieve L2 (genus-like) data in phyloseq format
pseq.L2 <- aggregate_taxa(pseq, level = "L2")
# Convert L2 to matrix format
genus.matrix.log10.simulated <- log10(otu_table(pseq.L2)@.Data)
```

```
## Error in eval(expr, envir, enclos): could not find function "otu_table"
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
# install.packages("googleVis") # install the library if needde

# See help(gvisMotionChart) for further details
mchart <- gvisMotionChart(df, idvar = "sample", timevar = "time")  
```

```
## Error in data[, vars.pos]: object of type 'closure' is not subsettable
```

```r
# Plot immediately (opens in browser, requires flash)
plot(mchart)  
```

```
## Error in plot(mchart): object 'mchart' not found
```

Save as html (needs javascript to open!). NOTE: html file viewing does not work locally - store the html file on server and view through internet:


```r
print(mchart, file="~/MotionChart.html")
```
