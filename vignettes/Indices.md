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
## Run 1 stress 0.253516 
## ... New best solution
## ... procrustes: rmse 0.007085436  max resid 0.09530197 
## Run 2 stress 0.2747136 
## Run 3 stress 0.2570603 
## Run 4 stress 0.2535947 
## ... procrustes: rmse 0.007484929  max resid 0.1099792 
## Run 5 stress 0.2540275 
## Run 6 stress 0.2536846 
## ... procrustes: rmse 0.006069665  max resid 0.06947634 
## Run 7 stress 0.253498 
## ... New best solution
## ... procrustes: rmse 0.006958935  max resid 0.1101667 
## Run 8 stress 0.2535004 
## ... procrustes: rmse 0.001213696  max resid 0.02143185 
## Run 9 stress 0.254337 
## Run 10 stress 0.2541912 
## Run 11 stress 0.2537164 
## ... procrustes: rmse 0.006747147  max resid 0.1105163 
## Run 12 stress 0.2536039 
## ... procrustes: rmse 0.003975626  max resid 0.07051828 
## Run 13 stress 0.2549219 
## Run 14 stress 0.2559363 
## Run 15 stress 0.253836 
## ... procrustes: rmse 0.01044807  max resid 0.1218678 
## Run 16 stress 0.2549497 
## Run 17 stress 0.2538858 
## ... procrustes: rmse 0.01152555  max resid 0.1217394 
## Run 18 stress 0.2540276 
## Run 19 stress 0.2552395 
## Run 20 stress 0.2738907
```

```
## Error in eval(expr, envir, enclos): object 'age' not found
```
