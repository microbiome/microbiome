## Intestinal microbiota diversity in 1006 western adults

Let us investigate an example data set from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). This contains large-scale profiling of 130 genus-like taxa across 1006 normal western adults. The data set is available from the [Data Dryad](http://doi.org/10.5061/dryad.pk75d) repository.


### Download HITChip Atlas data

[Load the HITChip Atlas microbiome profiling data in R](Data.md)


```r
# Download the required R packages and then the HITChip Atlas data set
library("rdryad")
library("microbiome")
d <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```


### Diversity 

# Estimate microbial diversity with different diversity measures


```r
library(phyloseq)
div <- estimate_richness(d, measures = c("Chao1", "Shannon"))
```

Plot richness vs. BMI with phyloseq tools (assuming you have the bmi_group field available in the sample metadata)


```r
plot_richness(d, x = "bmi_group", measures = c("Chao1", "Shannon"))
```

![plot of chunk div-example2](figure/div-example2-1.png) 


### Compare with known background factors

Diversity vs. obesity


```r
# Add diversity into sample data
sample_data(d)$diversity <- estimate_richness(d, measures = c("Shannon"))$Shannon

# Remove samples with no BMI information
dsub <- subset_samples(d, !is.na(bmi_group))

# Visualize
df <- harmonize_fields(sample_data(dsub))
p <- ggplot(df)
p <- p + geom_boxplot(aes(x = bmi_group, y = diversity))
p <- p + ggtitle("Microbiota diversity vs. obesity")
print(p)
```

![plot of chunk diversitywithmetadata](figure/diversitywithmetadata-1.png) 


Diversity vs. age with smoothed confidence intervals:


```r
library(microbiome)
library(sorvi)

# Pick the subset of RBB-preprocessed samples from time point 0
dsub <- subset_samples(d, time == 0 & DNA_extraction_method == "r")

# Collect variables into a data frame
df <- harmonize_fields(sample_data(dsub))

# Visualize
p <- sorvi::regression_plot(diversity~age, df)
print(p)
```

![plot of chunk atlas-example3](figure/atlas-example3-1.png) 


## Further resources

For further examples, see [microbiome tutorial](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)
