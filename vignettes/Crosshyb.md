## Visualizing cross-hybridization

To visualize cross-hybridization between selected taxa on HITChip (or
other chips), use the following.

By default the groups with no cross-hyb are filtered out for clarity. Rows and columns are ordered by hierarchical clustering. The cross-hyb is shown in percentages, rounded as indicated by the rounding argument. The percentage indicates which fraction of probes for a row taxon overlaps with probes of a column taxon. This is not symmetric if the row and col taxon have a different total number of probes. For details, see help(PlotCrosshyb).


```r
library(microbiome, quietly = TRUE)

# Check cross-hyb between all L1 groups
res <- PlotCrosshyb(tax.level = "L1", rounding = 1, show.plot = FALSE)
```

```
## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.filtered.tab
## No cross-hybriziation at L1 level
```

```r
# Pick the crosshyb table and figure
crosshyb.table <- res$data
p <- res$plot

# Plot the figure    
print(p)
```

```
## NULL
```

```r
# Organize the Crosshyb table
s <- subset(res$data, crosshyb > 0)
```

```
## Error in subset.default(res$data, crosshyb > 0): object 'crosshyb' not found
```

```r
s <- s[rev(order(s$crosshyb)),]
```

```
## Error in eval(expr, envir, enclos): object 's' not found
```

```r
head(s)
```

```
## Error in head(s): error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 's' not found
```

```r
#write.table(s, file = "crosshyb.tab", sep = "\t", quote = F, row.names = F)
```


### Further examples

Investigate species-species cross-hybridization within the Dialister L2 group


```r
# Pick the phylogeny which was used to summarize probes to species level
phylogeny.info <- GetPhylogeny("HITChip", "filtered") 
```

```
## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.filtered.tab
```

```r
# Select species belonging to Dialister L2 group
mytaxa <- levelmap("Dialister", level.from = "L2", level.to = "species", phylogeny.info = phylogeny.info)[[1]]

# Check cross-hyb between Dialister species
res <- PlotCrosshyb(tax.level = "species", selected.taxa = mytaxa, rounding = 0, phylogeny.info = phylogeny.info)
```

```
## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.filtered.tab
## The "ward" method has been renamed to "ward.D"; note new "ward.D2"
## The "ward" method has been renamed to "ward.D"; note new "ward.D2"
```

![plot of chunk chyb2](figure/chyb2-1.png) 

```r
# Check the cross-hyb data as well
head(res$data)
```

```
##                                   Taxon1                 Taxon2  crosshyb
## 1                      Dialister invisus      Dialister invisus   0.00000
## 2                 Dialister pneumosintes      Dialister invisus  20.00000
## 3 Uncultured bacterium clone Eldhufec089      Dialister invisus 100.00000
## 4 Uncultured bacterium clone Eldhufec093      Dialister invisus 100.00000
## 5              uncultured bacterium MG10      Dialister invisus   0.00000
## 6                      Dialister invisus Dialister pneumosintes  33.33333
```

