## HITChip taxonomy

Check the overall HITChip taxonomy:


```r
require(microbiome)
data("hitchip.taxonomy")
tax.table <- hitchip.taxonomy$full
kable(head(tax.table))
```



|    |L1             |L2               |species                |specimen               |oligoID  |L0             |
|:---|:--------------|:----------------|:----------------------|:----------------------|:--------|:--------------|
|sp1 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1134 |Actinobacteria |
|sp2 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1158 |Actinobacteria |
|sp3 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1194 |Actinobacteria |
|sp4 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1589 |Actinobacteria |
|sp5 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 1590 |Actinobacteria |
|sp6 |Actinobacteria |Actinomycetaceae |Actinomyces naeslundii |Actinomyces naeslundii |HIT 5644 |Actinobacteria |

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


