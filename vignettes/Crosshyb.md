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
## $colors
##  [1] "#0000FF" "#0C0CFF" "#1919FF" "#2626FF" "#3333FF" "#3F3FFF" "#4C4CFF"
##  [8] "#5959FF" "#6666FF" "#7272FF" "#7F7FFF" "#8C8CFF" "#9999FF" "#A5A5FF"
## [15] "#B2B2FF" "#BFBFFF" "#CCCCFF" "#D8D8FF" "#E5E5FF" "#F2F2FF" "#FFFFFF"
## [22] "#FFF2F2" "#FFE5E5" "#FFD8D8" "#FFCBCB" "#FFBFBF" "#FFB2B2" "#FFA5A5"
## [29] "#FF9898" "#FF8C8C" "#FF7F7F" "#FF7272" "#FF6565" "#FF5959" "#FF4C4C"
## [36] "#FF3F3F" "#FF3232" "#FF2626" "#FF1919" "#FF0C0C" "#FF0000"
## 
## $breaks
##  [1] -1000001.90       -1.95       -1.85       -1.75       -1.65
##  [6]       -1.55       -1.45       -1.35       -1.25       -1.15
## [11]       -1.05       -0.95       -0.85       -0.75       -0.65
## [16]       -0.55       -0.45       -0.35       -0.25       -0.15
## [21]       -0.05        0.05        0.15        0.25        0.35
## [26]        0.45        0.55        0.65        0.75        0.85
## [31]        0.95        1.05        1.15        1.25        1.35
## [36]        1.45        1.55        1.65        1.75        1.85
## [41]        1.95  1000001.90
## 
## $palette.function
## function (n) 
## {
##     x <- ramp(seq.int(0, 1, length.out = n))
##     if (ncol(x) == 4L) 
##         rgb(x[, 1L], x[, 2L], x[, 3L], x[, 4L], maxColorValue = 255)
##     else rgb(x[, 1L], x[, 2L], x[, 3L], maxColorValue = 255)
## }
## <bytecode: 0x163e1a40>
## <environment: 0x38ce8c90>
```

```r
# Organize the Crosshyb table
suppressMessages(library(dplyr))
s <- filter(res$data, crosshyb > 0)
```

```
## Error in UseMethod("filter_"): no applicable method for 'filter_' applied to an object of class "phyloseq"
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
## Error in x[seq_len(n)]: object of type 'S4' is not subsettable
```

