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
library(dplyr)

# Pick the phylogeny which was used to summarize probes to species level
tax.table <- GetPhylogeny("HITChip", "full")
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
#data("hitchip.taxonomy")
#tax.table <- as.data.frame(hitchip.taxonomy$filtered@.Data)

# Check cross-hyb between all L2 groups
res <- PlotCrosshyb(tax.level = "L2", rounding = 1, show.plot = FALSE, tax.table = tax.table)
```

```
## Error in tax_table(as(x, "matrix")[i, j, drop = FALSE]): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error in as(x, "matrix")[i, j, drop = FALSE] : subscript out of bounds
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
# Select species belonging to Dialister L2 group
mytaxa <- levelmap("Dialister", from = "L2", to = "species", tax.table)[[1]]
```

```
## Error in tax_table(as(x, "matrix")[i, j, drop = FALSE]): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error in as(x, "matrix")[i, j, drop = FALSE] : subscript out of bounds
```

```r
# Check cross-hyb between Dialister species
res <- PlotCrosshyb(tax.level = "species", selected.taxa = mytaxa, rounding = 0, tax.table = tax.table)
```

```
## Error in tax_table(as(x, "matrix")[i, j, drop = FALSE]): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error in as(x, "matrix")[i, j, drop = FALSE] : subscript out of bounds
```

```r
# Check the cross-hyb data as well
library(knitr)
kable(head(res$data))
```

```
## Error in kable_markdown(x = structure(character(0), .Dim = c(0L, 0L), .Dimnames = list(: the table must have a header (column names)
```

