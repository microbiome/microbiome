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
## Downloading data set from Lahti et al. Nat. Comm. 5:4344, 2014 from 
##   		       Data Dryad: http://doi.org/10.5061/dryad.pk75d
```

```
## Error in if (format == "phyloseq") {: argument is of length zero
```


### Diversity 

### Estimating microbial diversity with different diversity measures


```r
library(phyloseq)
div <- estimate_diversity(pseq, measures = c("Observed", "Shannon", "Simpson"))
kable(head(div))
```



|         | Observed|   Shannon|   Simpson|
|:--------|--------:|---------:|---------:|
|Sample.1 |       16| 1.5480650| 0.6362939|
|Sample.2 |       16| 0.6408117| 0.2881798|
|Sample.3 |       16| 0.6171970| 0.2854109|
|Sample.4 |       16| 1.7280468| 0.6839380|
|Sample.5 |       16| 0.6235029| 0.2930380|
|Sample.6 |       16| 0.4597322| 0.2264956|


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
