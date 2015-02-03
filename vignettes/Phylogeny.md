### Phylogeny

Map phylotypes between hierarchy levels:


```r
require(microbiome)
phylogeny.info <- GetPhylogeny("HITChip")
m <- levelmap(phylotypes = NULL, 
              level.from = "species", 
	      level.to = "L2", 
	      phylogeny.info = phylogeny.info)
```

