## HITChip taxonomy

Check the overall HITChip taxonomy:


```r
require(microbiome)
tax.table <- GetPhylogeny("HITChip")
head(tax.table)
```

```
## Taxonomy Table:     [6 taxa by 5 taxonomic ranks]:
##     L1               L2                 species                 
## sp1 "Actinobacteria" "Actinomycetaceae" "Actinomyces naeslundii"
## sp2 "Actinobacteria" "Actinomycetaceae" "Actinomyces naeslundii"
## sp3 "Actinobacteria" "Actinomycetaceae" "Actinomyces naeslundii"
## sp4 "Actinobacteria" "Actinomycetaceae" "Actinomyces naeslundii"
## sp5 "Actinobacteria" "Actinomycetaceae" "Actinomyces naeslundii"
## sp6 "Actinobacteria" "Actinomycetaceae" "Actinomyces naeslundii"
##     specimen                 oligoID   
## sp1 "Actinomyces naeslundii" "HIT 1134"
## sp2 "Actinomyces naeslundii" "HIT 1158"
## sp3 "Actinomyces naeslundii" "HIT 1194"
## sp4 "Actinomyces naeslundii" "HIT 1589"
## sp5 "Actinomyces naeslundii" "HIT 1590"
## sp6 "Actinomyces naeslundii" "HIT 5644"
```

Conversion between taxonomic levels:


```r
m <- levelmap(phylotypes = c("Akkermansia", "Bacteroides fragilis et rel."), 
              from = "L2", to = "L1", tax.table)

# Another example
data(GlobalPatterns)
taxtable <- tax_table(GlobalPatterns)
levelmap("Crenarchaeota", 'Phylum', 'Kingdom', taxtable)
```

```
## Taxonomy Table:     [1 taxa by 1 taxonomic ranks]:
##        Kingdom  
## 549322 "Archaea"
```


