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
|sample.1 |        5| 1.332441| 0.6886460|
|sample.2 |        5| 1.474480| 0.7391192|
|sample.3 |        5| 1.523999| 0.7624925|
|sample.4 |        5| 1.560578| 0.7804036|
|sample.5 |        5| 1.544700| 0.7721668|
|sample.6 |        5| 1.542622| 0.7745688|


### Diversity vs. obesity


```r
p <- plot_diversity(pseq, x = "bmi_group", measures = c("Observed", "Shannon", "Simpson"), det.th = 250)
```

```
## Error in `$<-.data.frame`(`*tmp*`, "horiz", value = structure(integer(0), .Label = character(0), class = "factor")): replacement has 0 rows, data has 44
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
## Error in eval(expr, envir, enclos): object 'DNA_extraction_method' not found
```

```r
# Visualize
library(sorvi)
p <- sorvi::regression_plot(diversity~age, sample_data(pseq))
```

```
## Error in eval(expr, envir, enclos): incorrect size (1), expecting : 44
```

```r
print(p)
```

![plot of chunk atlas-example3](figure/atlas-example3-1.png)


## Further resources

For further examples, see [microbiome tutorial](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)
