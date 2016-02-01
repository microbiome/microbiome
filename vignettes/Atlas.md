## Intestinal microbiota diversity in 1006 western adults

Let us investigate an example data set from [Lahti et al. Nat. Comm. 5:4344, 2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html). This contains microbiota profiling of 130 genus-like taxa across 1006 normal western adults from [Data Dryad](http://doi.org/10.5061/dryad.pk75d).


### Download HITChip Atlas data

[Load the HITChip Atlas microbiome profiling data in R](Data.md)


```r
# Download the required R packages and then the HITChip Atlas data set
library("rdryad")
library("microbiome")
pseq <- download_microbiome("atlas1006")
```

```
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```
## Error in curl::curl_fetch_memory(url, handle = handle): Timeout was reached
```


### Diversity 

### Estimating microbial diversity with different diversity measures


```r
library(phyloseq)
div <- estimate_diversity(pseq, measures = c("Observed", "Shannon", "Simpson"))
kable(head(div))
```



|         | Observed|  Shannon|   Simpson|
|:--------|--------:|--------:|---------:|
|Sample.1 |      130| 2.942723| 0.8677771|
|Sample.2 |      130| 2.824184| 0.8766237|
|Sample.3 |      130| 2.409584| 0.7670464|
|Sample.4 |      130| 2.994672| 0.8740136|
|Sample.5 |      130| 2.108225| 0.6847733|
|Sample.6 |      130| 2.073329| 0.6616150|


### Diversity vs. obesity


```r
p <- plot_diversity(pseq, x = "bmi_group", measures = c("Observed", "Shannon", "Simpson"), det.th = 250)
```

```
## Error in access(object, "otu_table", errorIfNULL): otu_table slot is empty.
```

```r
print(p)
```

![plot of chunk div-example2](figure/div-example2-1.png)


### Diversity vs. age


```r
# Pick the subset of RBB-preprocessed samples from time point 0
pseq <- subset_samples(pseq, time == 0 & DNA_extraction_method == "r")
```

```
## Error in time == 0: comparison (1) is possible only for atomic and list types
```

```r
# Visualize
library(sorvi)
p <- sorvi::regression_plot(diversity~age, sample_data(pseq))
```

```
## Error in eval(expr, envir, enclos): incorrect size (1), expecting : 222
```

```r
print(p)
```

![plot of chunk atlas-example3](figure/atlas-example3-1.png)


## Further resources

For further examples, see [microbiome tutorial](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)
