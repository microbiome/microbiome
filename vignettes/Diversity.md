### Diversity estimation


```r
# Determine detection threshold as the top 80 percent quantile
# of the data
det.th <- quantile(oligo.data, 0.8)

# Visualize the detection threshold (at log10 scale for clarity)
plot(density(log10(oligo.data))); abline(v = log10(det.th), main = "Detection threshold", xlab = "Abundance (Log10)", ylab = "Frequency")
```

![plot of chunk diversity-example](figure/diversity-example-1.png) 

```r
# Calculate richness. 
# This indicates how many oligos are present in each sample
# (exceed the detection threshold)
ri <- colSums(oligo.data > det.th)

# Diversity using the vegan package
# NOTE: data needs to be in absolute scale, not logarithmic
di <- vegan::diversity(t(oligo.data), index = "shannon")

# Pielou's evenness is S/ln(R) w.r.t. given detection threshold
# NOTE: here we use detection threshold for diversity as well because
# the exact same data has to be used for diversity and richness calculation,
# and for richness calculation the detection threshold needs to be set anyway
# Diversity can be as such calculated also without threshold (see above)
# but that gives somewhat different result.
oligo.data2 <- oligo.data - det.th # NOTE: absolute (not log) scale data
S <- vegan::diversity(t(oligo.data2), index = "shannon")
R <- colSums(oligo.data2 > 0)
ev <- S/log(R)

# Combine all into a single table
divtab <- cbind(richness = ri, evenness = ev, diversity = di)
head(divtab)
```

```
##          richness  evenness diversity
## Sample.1      918 0.9076456  6.015876
## Sample.2      728 0.9243413  5.759365
## Sample.3      990 0.8960132  6.137286
## Sample.4      698 0.8896066  5.630066
## Sample.5      661 0.8952997  5.595988
## Sample.6      334 1.0606226  5.327884
```


Calculate diversities within each phylum-level group


```r
L1.groups <- unique(phylogeny.info.filtered$L1) # List all L1 groups
L1.diversities <- NULL # initialize
for (phylum in L1.groups) {
  # List all probes for the given phylum
  probes <- levelmap(phylum, "L1", "oligoID", phylogeny.info = phylogeny.info.filtered)[[phylum]]
  # Check diversity within this phylum
  di <- vegan::diversity(t(oligo.data[probes,]), index = "shannon")
  # Add to the table
  L1.diversities <- cbind(L1.diversities, di)
}
# Name the columns
colnames(L1.diversities) <- L1.groups
```


### Diversity boxplot

Diversity boxplot for selected samples:


```r
group <- metadata$group
names(group) <- metadata$sampleID
my.samples <- names(group)
boxplot(di[my.samples]  ~ group[my.samples], las = 1)
```

![plot of chunk diversity-example2](figure/diversity-example2-1.png) 

