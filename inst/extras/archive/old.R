### Comparing of two or more groups with a parametric test (linear model; ANOVA)

Note: check ANOVA modeling assumptions before testing. 

```{r comparisons2}
# Load example data from a 
# [diet swap study](http://dx.doi.org/10.1038/ncomms7342)
data("dietswap")
pseq <- dietswap

# Convert to relative abundances
pseq <- transform_phyloseq(pseq, "relative.abundance")

# 1-way ANOVA p-values for the multi-group comparison across time groups
source(system.file("extdata/check_anova.R", package = "microbiome"))
anova.results <- check_anova(pseq, "group", p.adjust.method = "BH")
kable(head(anova.results))
```

### Wilcoxon test (two-group comparisons)

If the data remarkably violates Gaussian assumptions, try
non-parametric test. Wilcoxon is one option for two group
comparison. Here we compare males and females in the example data. The
check_wilcoxon function is experimental and not included in the
microbiome R package but is available as a supplementary file:

```{r comparisons-exampless}
source(system.file("extdata/check_wilcoxon.R", package = "microbiome"))
pval <- check_wilcoxon(pseq, "sex")
```

### Comparing limma and t-test

Order the taxa with t-test for comparison and validation purposes. The
differences are small in this simulated example, but [can be
considerable in real
data](http://www.plosone.org/article/info:doi/10.1371/journal.pone.0012336).

```{r limma-compairson, warning=FALSE, fig.path = "figure/"}
# Compare the two groups with t-test
library(dplyr)
pvalues.ttest <- c()
male.samples <- suppressWarnings(dplyr::filter(meta, gender == "male")$sample)
female.samples <- suppressWarnings(dplyr::filter(meta, gender == "female")$sample)
for (tax in rownames(otu)) {
  pvalues.ttest[[tax]] <- t.test(otu[tax, male.samples], otu[tax, female.samples])$p.value
}
# Multiple testing correction
pvalues.ttest <- p.adjust(pvalues.ttest, method = "fdr")

# Compare p-values between limma and t-test
taxa <- rownames(otu)
plot(pvalues.ttest[taxa], pvalues.limma[taxa])
abline(0,1,lty = 2)
```



#' read_profiling
#' 
#' Read run.profiling.script output into R
#'
#' @param data.dir Profiling script output directory for reading the data. 
#'                 If not given, GUI will ask to specify the file and 
#'             	   overruns the possible level / method arguments in the 
#'             	   function call.
#' @param method Select the preprocessing method that you like to check
#' 
#' @return data matrix (phylo x samples)
#'
#' @export
#' @examples 
#'   # data.dir <- system.file("extdata", package = "microbiome")
#'   # dat <- read_profiling(data.dir)
#'
#' @references See citation('microbiome')
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
read_profiling <- function(data.dir, method) {

  message(paste("Reading Chip data from", data.dir))
  results <- list()

  # Read probe-level data
  f <- paste(data.dir, "/oligoprofile.tab", sep = "")
  tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
  colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
  results[["probedata"]] <- tab

  # Read abundance tables
  for (s in c("L1", "L2", "species")) {
    f <- paste(data.dir, "/species-", method, ".tab", sep = "")
    tab <- read.csv(f, header = TRUE, sep = "\t", row.names = 1, as.is = TRUE)
    colnames(tab) <- unlist(strsplit(readLines(f, 1), "\t"))[-1]
    results[[s]] <- tab
  }
  
  # Read taxonomy table
  f <- paste(data.dir, "/taxonomy.tab", sep = "")
  if (!file.exists(f)) {
    f <- paste(data.dir, "/phylogeny.filtered.tab", sep = "")  
  }
  taxonomy <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  results[["taxonomy"]] <- tab
    
  # Read unfiltered taxonomy table
  f <- paste(data.dir, "/taxonomy.full.tab", sep = "")
  if (!file.exists(f)) {
    f <- paste(data.dir, "/phylogeny.full.tab", sep = "")  
  }  
  taxonomy.full <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
  results[["taxonomy.full"]] <- tab

  # Read sample metadata      
  f <- paste(data.dir, "/meta.tab", sep = "")
  if (file.exists(f)) {
    tab <- read.csv(f, header = TRUE, sep = "\t", as.is = TRUE)
    rownames(tab) <- tab$sample
    meta <- tab
    results[["meta"]] <- tab
  }

  results
    
} 


