

## Group-wise comparisons

Let us first read simulated example data.  See
[read.profiling](reading) on how to use your own data. 


```r
# Load the package
library(microbiome)  

# Load data (here simulated example data) 
data.directory <- system.file("extdata", package = "microbiome")
hitchip.data <- read.profiling(level = "L2", method = "frpa", data.dir = data.directory, log10 = TRUE)  

# Read simulated example metadata
library(gdata)
metadata.file <- paste(data.directory, "/metadata.xls", sep = "")
metadata <- read.xls(metadata.file, as.is = TRUE)
rownames(metadata) <- metadata$sampleID

# Ensure that the samples are in the same order in 
# HITChip data and metadata
metadata <- metadata[colnames(hitchip.data),]
head(metadata)
```

```
##          sampleID time      age      bmi   subjectID   group gender
## Sample.1 Sample.1    1 85.84339 32.28638 subjectID-1 group-3      M
## Sample.2 Sample.2    2 94.62583 23.11057 subjectID-2 group-1      M
## Sample.3 Sample.3    3 99.95209 36.18861 subjectID-3 group-1      F
## Sample.4 Sample.4    4 20.30884 31.63783 subjectID-4 group-3      M
## Sample.5 Sample.5    1 23.24289 28.46594 subjectID-1 group-2      F
## Sample.6 Sample.6    2 58.36422 28.22124 subjectID-2 group-1      F
##               diet  X
## Sample.1 Beverages NA
## Sample.2 Beverages NA
## Sample.3   Carrots NA
## Sample.4   Carrots NA
## Sample.5  Apricots NA
## Sample.6  Apricots NA
```

### Comparing of two or more groups with a parametric test (linear model; ANOVA)

NOTE: the HITChip data (hitchip.data) must be at log10 scale!


```r
# Define sample groups based on a selected metadata column
group <- metadata$diet
names(group) <- metadata$sampleID

# Calculate p-values for multi-group comparison
pvals <- check.anova(hitchip.data, group, p.adjust.method = "BH", sort = TRUE)

# Write the results to a table
write.table(pvals, file = "myfile.tab", quote = F)
```


### Wilcoxon test (only for two-group comparison)

If the data remarkably violates Gaussian assumptions, you need to use
non-parametric test, Wilcoxon is one option for two group
comparison. Here we compare males and females in simulated toy
data. The sampleIDs in G1 and G2 must match with input matrix
columns. Compare two sample groups and return adjusted p-values for
each phylotype. This test uses unpaired Wilcoxon test with BH
correction.


```r
# Define your sample groups
males <- subset(metadata, gender == "M")$sampleID
females <- subset(metadata, gender == "F")$sampleID

# Calculate Wilcoxon test with BH p-value correction for genus.matrix
pval <- check.wilcoxon(hitchip.data, G1 = males, G2 = females, p.adjust.method = "BH", sort = TRUE)  

# Investigate the top findings
head(pval)
```


