
To run these examples, you need to download the [HITChip Atlas data set](Data.md)

### Richness 


```r
library(microbiome)

# Get some HITChip data in phyloseq format
pseq <- download_microbiome("atlas1006")

# Pick the OTU data
# (note the zero point has been moved to the detection threshold;
#  typically signal 1.8 at HITChip log10 scale)
otu <- otu_table(pseq)@.Data

# Determine detection threshold at the 0.15 quantile
# Bacteria that exceed this threshold are considered present
# otherwise absent
det.th <- quantile(otu, 0.15)

# Visualize the detection threshold (at log10 scale for clarity)
plot(density(log10(otu)), main = "Detection threshold", xlab = "Abundance (Log10)", ylab = "Frequency")
abline(v = log10(det.th))
```

![plot of chunk rich-example](figure/rich-example-1.png) 

```r
# Calculate richness.
# This indicates how many oligos are present in each sample
# (exceed the detection threshold)
ri <- rowSums(otu > det.th)
hist(ri, main = "Richness")
```

![plot of chunk rich-example](figure/rich-example-2.png) 


Highlight specific groups:


```r
library(ggplot2)
data.dietswap <- download_microbiome("dietswap")
p <- plot_richness(data.dietswap, x = "gender", color = "group", measures = c("Shannon", "Simpson")) 
p <- p + geom_boxplot()
```


### Diversity 



```r
di <- estimate_richness(pseq, measures = c("Shannon"))
hist(di$Shannon, main = "Diversity")
```

![plot of chunk div-example](figure/div-example-1.png) 

Plot richness vs. BMI with phyloseq tools (assuming you have the bmi_group field available in the sample metadata)


```r
p <- plot_richness(pseq, x = "bmi_group", measures = c("Chao1", "Shannon"))
p <- p + geom_boxplot()
print(p)
```

![plot of chunk div-example2](figure/div-example2-1.png) 

### Compare with known background factors

Diversity vs. obesity group (factor):


```r
p <- plot_diversity(pseq, x = "bmi_group", measures = c("Shannon", "Simpson"))
print(p)
```

![plot of chunk diversitywithmetadata](figure/diversitywithmetadata-1.png) 

Diversity vs. age (continuous):


```r
p <- plot_diversity(pseq, x = "age", measures = "Shannon")
print(p)
```

![plot of chunk diversitywithmetadata2](figure/diversitywithmetadata2-1.png) 


Diversity vs. age with smoothed confidence intervals - manual version:


```r
library(microbiome)
library(sorvi)
library(dplyr)

# Add diversity into sample metadata
sample_data(pseq)$diversity <- estimate_richness(pseq, measures = c("Shannon"))$Shannon

# Select a subset of samples
pseq0 <- subset_samples(pseq, time == 0 & DNA_extraction_method == "r")

# Visualize
df <- sample_data(pseq0)
p <- sorvi::regression_plot(diversity ~ age, df)
print(p)
```

![plot of chunk diversity-example13](figure/diversity-example13-1.png) 

