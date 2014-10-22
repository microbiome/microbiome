# Normalization

Normalization is used to correct for technical differences in sample
intensities. Min/max is the recommended default methods; quantile normalization should be used only when it can be assumed that the overall distribution of phylotypes (the phylotype identities can vary, though) within each sample is approximately similar. These standard options correspond to different background assumptions:

* **Min/max** The default option. Aims to normalize the samples
to comparable scales without affecting the phylotype distribution
within a sample. Particularly useful when there are considerable
systematic differences between samples, for instance due to
antiobiotic treatment. Sets the minimum and maximum values of each
sample to the same value by shifting and scaling the samples
accordingly. These values are selected as the 0.5% 99.5% quantiles
of the data to improve robustness.

* **Quantile** The quantile normalization forces the same phylotype distribution on all samples, which is calculated based on the average over all samples. Quantile normalization may be more efficient in removing technical biases from the data than min/max but it assumes that the overall phylotype distribution is approximately same for the different samples (ie. there are no antiobiotic treatments, and the samples/subjects are also otherwise similar). 


* **None** skip normalization 

