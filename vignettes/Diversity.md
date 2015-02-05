### Diversity estimation

Download the [HITChip Atlas data set](Data.md)


```r
# Use the HITChip Atlas data
mydata <- atlas$microbes

# Determine detection threshold as the top-80 quantile
det.th <- quantile(mydata, 0.8)

# Visualize the detection threshold (at log10 scale for clarity)
plot(density(log10(mydata)), main = "Detection threshold", xlab = "Abundance (Log10)", ylab = "Frequency")
abline(v = log10(det.th))
```

![plot of chunk diversity-example](figure/diversity-example-1.png) 

```r
# Calculate richness. 
# This indicates how many oligos are present in each sample
# (exceed the detection threshold)
ri <- rowSums(mydata > det.th)
hist(ri, main = "Richness")
```

![plot of chunk diversity-example](figure/diversity-example-2.png) 

```r
# Diversity using the vegan package
# NOTE: data needs to be in absolute scale, not logarithmic
di <- vegan::diversity(mydata, index = "shannon")
hist(di, main = "Diversity")
```

![plot of chunk diversity-example](figure/diversity-example-3.png) 

```r
# Pielou's evenness is S/ln(R) w.r.t. given detection threshold
# NOTE: here we use detection threshold for diversity as well because
# the exact same data has to be used for diversity and richness calculation,
# and for richness calculation the detection threshold needs to be set anyway
# Diversity can be as such calculated also without threshold (see above)
# but that gives somewhat different result.
mydata2 <- mydata - det.th # NOTE: absolute (not log) scale data
S <- vegan::diversity(mydata2, index = "shannon")
R <- rowSums(mydata2 > 0)
ev <- S/log(R)

# Combine all into a single table
divtab <- cbind(richness = ri, evenness = ev, diversity = di)
head(divtab)
```

```
##          richness  evenness diversity
## Sample-1       25  7.758274  3.249386
## Sample-2       34  1.544166  3.438738
## Sample-3       13  7.186927  2.935857
## Sample-4       25  1.161892  3.103601
## Sample-5       27  1.216306  3.119424
## Sample-6       18 13.354305  3.007597
```


### Compare diversity with known background factors

Start by loading the [HITChip Atlas data set and metadata](Data.Rmd)


```r
par(mar = c(6, 4, 3, 1), mfrow = c(1, 2))
d <- atlas$microbes
bmi <- atlas$meta$BMI_group
age <- atlas$meta$Age
di <- vegan::diversity(d)
boxplot(di ~ bmi, las = 2, main = "Microbiota diversity vs. obesity")
plot(age, di, main = "Microbiota diversity vs. Age", ylab = "Diversity", xlab = "Age (years)")
```

![plot of chunk diversitywithmetadata](figure/diversitywithmetadata-1.png) 


Plot subject age versus phylotype abundance with smoothed confidence intervals:


```r
library(microbiome)
library(sorvi)
library(dplyr)

# Just pick RBB-preprocessed samples from time point 0
rbb.samples <- filter(atlas$meta, Time == 0 & DNA_extraction_method == "r")$SampleID

# Collect variables into a data frame
df <- data.frame(Age = atlas$meta[rbb.samples, "Age"], Diversity = di[rbb.samples])

# Regression plot
p <- sorvi::regression_plot(Diversity~Age, df, shade = TRUE, mweight = TRUE, verbose = FALSE)
print(p)
```

![plot of chunk visu-example3](figure/visu-example3-1.png) 

