## HITChip Phylogeny

Map HITChip phylogenetic microarray phylotypes between taxonomic hierarchy levels.  


Check the overall phylogeny table:


```r
require(microbiome)
phylogeny.info <- GetPhylogeny("HITChip")
head(phylogeny.info)
```

```
##               L1               L2                species
## 1 Actinobacteria Actinomycetaceae Actinomyces naeslundii
## 2 Actinobacteria Actinomycetaceae Actinomyces naeslundii
## 3 Actinobacteria Actinomycetaceae Actinomyces naeslundii
## 4 Actinobacteria Actinomycetaceae Actinomyces naeslundii
## 5 Actinobacteria Actinomycetaceae Actinomyces naeslundii
## 6 Actinobacteria Actinomycetaceae Actinomyces naeslundii
##                 specimen  oligoID
## 1 Actinomyces naeslundii HIT 1134
## 2 Actinomyces naeslundii HIT 1158
## 3 Actinomyces naeslundii HIT 1194
## 4 Actinomyces naeslundii HIT 1589
## 5 Actinomyces naeslundii HIT 1590
## 6 Actinomyces naeslundii HIT 5644
```

Convert taxa from one level to another level:


```r
m <- levelmap(phylotypes = c("Akkermansia", "Bacteroides fragilis et rel."), 
              level.from = "L2", 
	      level.to = "L1", 
	      phylogeny.info = phylogeny.info)

head(m)
```

```
## [1] Verrucomicrobia Bacteroidetes  
## Levels: Bacteroidetes Verrucomicrobia
```

