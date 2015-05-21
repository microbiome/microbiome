### Boxplots


```r
library(microbiome)
library(ggplot2)

# Microbiota profiling data. Read as: bacteria x samples matrix
# From https://peerj.com/articles/32/
pseq <- download_microbiome("peerj32")$physeq

p <- boxplot_abundance(pseq, x = "time", y = "Akkermansia", line = "subject", color = "gender")

print(p)
```

![plot of chunk boxplot-example](figure/boxplot-example-1.png) 
