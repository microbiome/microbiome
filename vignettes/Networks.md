
## Networks

Load example data:


```r
library(microbiome)
# Get some HITChip data in phyloseq format
pseq <- download_microbiome("dietswap")
```

```
## Error in curl::curl_fetch_memory(url, handle = handle): Timeout was reached
```

For more network examples, see [phyloseq tutorial](http://joey711.github.io/phyloseq/plot_network-examples)


```r
ig <- make_network(pseq, max.dist = 0.2)
plot_network(ig, pseq, color = "nationality", shape = "group", line_weight = 0.4, label = "sample")
```

```
## Error in plot_network(ig, pseq, color = "nationality", shape = "group", : The graph you provided, `g`, has too few vertices. 
##          Check your graph, or the output of `make_network` and try again.
```

Another example:


```r
plot_net(pseq, maxdist = 0.2, point_label = "group")
```

```
## Error in data.table(laymeth(g, ...), vertex = get.vertex.attribute(g, : column or argument 2 is NULL
```

