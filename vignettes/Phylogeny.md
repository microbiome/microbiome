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
## Taxonomy Table:     [6 taxa by 2 taxonomic ranks]:
##                              Phylum           
## Actinomycetaceae             "Actinobacteria" 
## Aerococcus                   "Bacilli"        
## Aeromonas                    "Proteobacteria" 
## Akkermansia                  "Verrucomicrobia"
## Alcaligenes faecalis et rel. "Proteobacteria" 
## Allistipes et rel.           "Bacteroidetes"  
##                              Genus                         
## Actinomycetaceae             "Actinomycetaceae"            
## Aerococcus                   "Aerococcus"                  
## Aeromonas                    "Aeromonas"                   
## Akkermansia                  "Akkermansia"                 
## Alcaligenes faecalis et rel. "Alcaligenes faecalis et rel."
## Allistipes et rel.           "Allistipes et rel."
```

Conversion between taxonomic levels:


```r
m <- levelmap(c("Akkermansia", "Bacteroides fragilis et rel."), 
              from = "L2", to = "L1", tax.table)
```

```
## Error in tax_table(as(x, "matrix")[i, j, drop = FALSE]): error in evaluating the argument 'object' in selecting a method for function 'tax_table': Error in as(x, "matrix")[i, j, drop = FALSE] : subscript out of bounds
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


