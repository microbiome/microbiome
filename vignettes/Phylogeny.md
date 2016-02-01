## HITChip taxonomy

Check the overall HITChip taxonomy:


```r
require(microbiome)
data("hitchip.taxonomy")
tax.table <- hitchip.taxonomy$full
kable(head(tax.table))
```



|    |L1             |L2               |species                |specimen               |oligoID  |
|:---|:--------------|:----------------|:----------------------|:----------------------|:--------|
|sp1 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1134 |
|sp2 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1158 |
|sp3 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1194 |
|sp4 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1589 |
|sp5 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1590 |
|sp6 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 5644 |

Conversion between taxonomic levels:


```r
m <- levelmap(c("Akkermansia", "Bacteroides fragilis et rel."), 
              from = "L2", to = "L1", tax.table)

# Another example
data(GlobalPatterns)
taxtable <- tax_table(GlobalPatterns)
levelmap("Crenarchaeota", 'Phylum', 'Kingdom', taxtable)
```

```
## [1] "Archaea"
```


