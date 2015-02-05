### Relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the
input data set needs to be in absolute scale (not logarithmic).


```r
rel <- relative.abundance(t(data))

# Rearrange the data for ggplot visualization tools
library(reshape)
dfm <- melt(rel)
colnames(dfm) <- c("Taxon", "SampleID", "RelativeAbundance")

# Provide barplot visualizations of relative abundances for some randomly selected samples
library(ggplot2)
dfmf <- filter(dfm, SampleID %in% c("Sample-1", "Sample-2", "Sample-3", "Sample-4", "Sample-5"))
p <- ggplot(dfmf, aes(x = SampleID, y = RelativeAbundance, fill = Taxon))
p <- p + geom_bar(position = "stack", stat = "identity")
print(p)
```

![plot of chunk diversity-example6](figure/diversity-example6-1.png) 

```r
# Also note that taking relative abundances likely changes the abundance histograms
theme_set(theme_bw(20))
p <- ggplot(filter(dfm, Taxon == "Prevotella.melaninogenica.et.rel."), aes(x = 100*RelativeAbundance))
p <- p + geom_density(fill = "darkgray")
p <- p + scale_x_log10()
p <- p + xlab("Relative Abundance (%)")
print(p)
```

![plot of chunk diversity-example6](figure/diversity-example6-2.png) 

