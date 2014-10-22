## Full preprocessing protocol for HITChip/MITChip/PITChip/ChickChip analysis

  1. Select folder for the output files 
  1. Select the studies
  1. Select samples (two or more)
  1. Use default parameters, or select to change the defaults: 
    * Select phylogeny if multiple options are available for your Chip (Default: 16S)
    * Select whether to remove (TRUE) or keep (FALSE) non-specific oligos (Default: FALSE)
    * Select [normalization method](normalization) (Default: minmax)  
  1. Wait until preprocessing completes (this can take 5-30 minutes)
  1. As the preprocessing completes, clustering results will be displayed for Euclidean metric and Pearson correlation for oligo-level data at absolute and log10 domains  
  1. Output files are in your selected output folder (from step 3), including:
  * Profiling log: details of preprocessing parameters in text file (\<date\>_\<time\>_profiling_log.txt) and in RData format (\<date\>_\<time\>_profiling_params.RData)
    * Phylotypes x samples data matrices. At absolute scale with alternative oligo summarization methods (sum, frpa) and taxonomic levels (species, L2, L1). For MITChip and PITChip, also L0 is provided. The data matrices are in the tab-separated \<level\>-\<method\>.tab files (for instance, L2-sum.tab). The 'sum' is the recommended default method. The 'nmf' is not available for species level. The 'ave' is not provided in output but can be retrieved afterwards with [preprocessing functions](preprocessing). For details, see [probe summarization](summarization). 
    * oligoprofile.tab (oligos vs. samples data matrix in the original scale; this file includes only those oligos that have passed the quality filtering steps)
    * phylogeny.tab (oligo-phylotype mappings and other oligo-specific information)
    * L1-oligoprofileClustering.png (oligos vs. samples heatmap; indicating the associated L1 groups). You can reproduce and modify the image further using the [PlotPhylochipHeatmap](oligoheatmap) function. If the sample size is high, plotting the image may fail and the image is empty.
    * hclust.png files: hierarchical clustering of the absolute and log10 oligo level data with euclidean and Pearson similarity measures
  1. If required, rerun the run.profiling.script to preprocess chip data with different parameters
  1. Copy the output directory (specified in step 3) on your own computer
  1. To read and analyse the preprocessed data further, see the [documentation page](analysis).


If needed, you can reproduce and modify the
L2-oligoprofileClustering.png (see
[PlotPhylochipHeatmap](oligoheatmap)) and hierarchical clustering (see
[calculate.hclust](clustering)) images generated during profiling
using the tools listed [here](visualization).




