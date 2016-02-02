## Microbiome-based indices

Various indices have been proposed for microbiome analyses. Here ready-made implementations for some of them.


### Maturity index

The microbiota maturity index has been adapted from [Subramanian, S. et al. (2014)](http://www.nature.com/nature/journal/v510/n7505/abs/nature13421.html) and [Dogra, S. et al. (2015)](http://mbio.asm.org/content/6/1/e02419-14.short), where microbiota maturity index has been shown to differentiate healthy children. In [Korpela et al. (2016)](www.nature.com/ncomms/2016/160126/ncomms10410/full/ncomms10410.html) this was calculated as the first principal coordinate from a PCoA (MDS) using only significantly age-associated genus-level taxa (all groups included). In this function NMDS is used instead. The maturity index is also adjusted for age.

This is an artificial example based on readily available adult data
set whereas the maturity index is typically calculated from and used
for babies/children:


```r
library(phyloseq)
data("atlas1006")
pseq <- atlas1006
pseq <- subset_samples(pseq, DNA_extraction_method == "r" & time == 0)
pseq <- transform_phyloseq(pseq, "relative.abundance")
maturity.index <- maturity(pseq)
```

```
## Error in eval(expr, envir, enclos): could not find function "maturity"
```
