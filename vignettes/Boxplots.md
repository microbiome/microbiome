### Boxplots


```r
# Load libraries
library(microbiome)
library(ggplot2)

# Probiotics intervention example data from https://peerj.com/articles/32/
data("peerj32")

# Abundance boxplot
p <- boxplot_abundance(peerj32$phyloseq, x = "time", y = "Akkermansia", line = "subject", color = "gender")
print(p)
```

![plot of chunk boxplot-example](figure/boxplot-example-1.png)
