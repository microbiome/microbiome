## HITChip taxonomy

Check the overall HITChip taxonomy:


```r
require(microbiome)
tax.table <- GetPhylogeny("HITChip")
```

```
## Error in file(file, "rt"): cannot open the connection
```

```r
head(tax.table)
```

```
## Error in head(tax.table): error in evaluating the argument 'x' in selecting a method for function 'head': Error: object 'tax.table' not found
```

Conversion between taxonomic levels:


```r
m <- levelmap(c("Akkermansia", "Bacteroides fragilis et rel."), 
              from = "L2", to = "L1", tax.table)
```

```
## Error in levelmap(c("Akkermansia", "Bacteroides fragilis et rel."), from = "L2", : object 'tax.table' not found
```

```r
# Another example
data(GlobalPatterns)
taxtable <- tax_table(GlobalPatterns)
levelmap("Crenarchaeota", 'Phylum', 'Kingdom', taxtable)
```

```
## [1] "Archaea"
```


