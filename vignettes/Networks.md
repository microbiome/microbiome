
## Networks

[Networks](http://joey711.github.io/phyloseq/plot_network-examples)


```r
library(microbiome)
data.dietswap <- download_microbiome("dietswap")

plot_net(physeq, maxdist = 0.45, point_label = "group")
```

```
## Error in distance(physeq, method, type): object 'physeq' not found
```

```r
ig <- make_network(physeq, max.dist = 0.45)
```

```
## Error in distance(physeq, method = distance, type = type, ...): object 'physeq' not found
```

```r
plot_network(ig, physeq, color = "gender", shape = "group", line_weight = 0.4, label = NULL)
```

```
## Error in match(x, table, nomatch = 0L): object 'ig' not found
```

