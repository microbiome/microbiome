
## Group-wise comparisons

Read example data from a [diet swap study](http://dx.doi.org/10.1038/ncomms7342):


```r
library(microbiome)
data("dietswap")
pseq <- dietswap
```

### Comparing of two or more groups with a parametric test (linear model; ANOVA)

Note that in practice it will be necessary to check ANOVA modeling assumptions before testing:


```r
# Convert to relative abundances
pseq <- transform_phyloseq(pseq, "relative.abundance")

# 1-way ANOVA p-values for the multi-group comparison across time groups
anova.results <- check_anova(pseq, "group", p.adjust.method = "BH")
kable(head(anova.results))
```



|   p.anova|   p.ED-DI|   p.HE-DI|   p.HE-ED|    ave.DI|    ave.ED|    ave.HE|
|---------:|---------:|---------:|---------:|---------:|---------:|---------:|
| 0.0010594| 0.9545821| 0.0001572| 0.0000392| 1.4442650| 1.3639240| 2.5845650|
| 0.0029303| 0.8320058| 0.0010983| 0.0001006| 0.5039235| 0.4327322| 0.9486485|
| 0.0336511| 0.9979732| 0.0038676| 0.0028199| 0.0376643| 0.0379573| 0.0219853|
| 0.0336511| 0.9646245| 0.0026158| 0.0054579| 0.1825854| 0.1921543| 0.3084474|
| 0.0706163| 0.0382707| 0.6660798| 0.0026381| 0.4153172| 0.7268781| 0.3066399|
| 0.0708257| 0.6389312| 0.0043033| 0.0533601| 0.4995820| 0.6407194| 1.0016500|


### Wilcoxon test (two-group comparisons)

If the data remarkably violates Gaussian assumptions use
non-parametric test. Wilcoxon is one option for two group
comparison. Here we compare males and females in the example data.


```r
pval <- check_wilcoxon(pseq, "sex")
```



