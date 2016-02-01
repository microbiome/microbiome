### TODO

 * Stability line plot ?

 * Tipping elements

 * Classification of the taxa into the distinct abundance types

 * Microbiota maturity index? Has previously been shown to differentiate between healthy children and those with health issues. Was calculated as the first principal coordinate from a PCoA using only significantly age-associated genus-level taxa in the control and early-life groups.
 
 * Negative binomial for sample comparisons

 * pairwise comparisons with subject random effect term

 * gower distance for sample similarities

 * SPIEC-EASE / SparCC network algorithms from the huge package -> Healthy vs other networks ?

 * Bagged and ordinary RDA for phyloseq objects with ggplot visualizations

 * Heatmap: indicate groups with a color bar

 * Homogeneity vs. sample size dependency. Normalize.

# Also PERMANOVA shows the communities are significantly different
# across groups
#adonis(t(dat.full) ~ group, data=meta.full, permutations=99)
