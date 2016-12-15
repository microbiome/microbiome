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
## Wisconsin double standardization
## Run 0 stress 0.2536304 
## Run 1 stress 0.2536459 
## ... procrustes: rmse 0.004656075  max resid 0.06960597 
## Run 2 stress 0.2534972 
## ... New best solution
## ... procrustes: rmse 0.0067322  max resid 0.07377297 
## Run 3 stress 0.254041 
## Run 4 stress 0.2539454 
## ... procrustes: rmse 0.01153968  max resid 0.123549 
## Run 5 stress 0.2540213 
## Run 6 stress 0.4189932 
## Run 7 stress 0.2536004 
## ... procrustes: rmse 0.003814867  max resid 0.07098161 
## Run 8 stress 0.254319 
## Run 9 stress 0.2539444 
## ... procrustes: rmse 0.009359277  max resid 0.1220958 
## Run 10 stress 0.2537099 
## ... procrustes: rmse 0.007066883  max resid 0.1097338 
## Run 11 stress 0.2535898 
## ... procrustes: rmse 0.002876847  max resid 0.05457971 
## Run 12 stress 0.2545886 
## Run 13 stress 0.4189923 
## Run 14 stress 0.2548786 
## Run 15 stress 0.2724054 
## Run 16 stress 0.2543242 
## Run 17 stress 0.2542938 
## Run 18 stress 0.2536374 
## ... procrustes: rmse 0.006842737  max resid 0.110375 
## Run 19 stress 0.4189799 
## Run 20 stress 0.2540882
```

```
## Error in eval(expr, envir, enclos): object 'age' not found
```
