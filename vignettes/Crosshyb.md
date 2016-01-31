## Visualizing cross-hybridization

To visualize cross-hybridization between selected taxa on HITChip (or
other chips), use the following scripts. By default the groups with no
cross-hyb are filtered out for clarity. Rows and columns are ordered
by hierarchical clustering. The cross-hyb is shown in percentages,
rounded as indicated by the rounding argument. The percentage
indicates which fraction of probes for a row taxon overlaps with
probes of a column taxon. This is not symmetric if the row and col
taxon have a different total number of probes. For details, see
help(PlotCrosshyb).


```r
library(microbiome, quietly = TRUE)

# Check cross-hyb between all L1 groups
res <- PlotCrosshyb(tax.level = "L2", rounding = 1, show.plot = FALSE)
```

```
## Error in match.fun(FUN): argument "FUN" is missing, with no default
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
suppressMessages(library(dplyr))
s <- filter(res$data, crosshyb > 0)
```

```
## Error in UseMethod("filter_"): no applicable method for 'filter_' applied to an object of class "NULL"
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


### Further examples

Investigate species-species cross-hybridization within the Dialister L2 group


```r
# Pick the phylogeny which was used to summarize probes to species level
tax.table <- GetPhylogeny("HITChip", "filtered") 

# Select species belonging to Dialister L2 group
mytaxa <- levelmap("Dialister", from = "L2", to = "species", tax.table)[[1]]

# Check cross-hyb between Dialister species
res <- PlotCrosshyb(tax.level = "species", selected.taxa = mytaxa, rounding = 0, tax.table)
```

```
## Error in match.fun(FUN): argument "FUN" is missing, with no default
```

```r
# Check the cross-hyb data as well
library(knitr)
kable(head(res$data))
```

```
## Error in kable_markdown(x = structure(character(0), .Dim = c(0L, 0L), .Dimnames = list(: the table must have a header (column names)
```

