# Overview of H/M/PITChip analysis

The [preprocessed data files](profiling) are stored into the output
directory which you specified during preprocessing. Copy the output
files on USB stick and analyse them further on any computer using R
with the [read.profiling](reading) function using the tools listed
below. 

Analysis of the H/M/PITChip data includes the following steps. 

 1. [Install the package](installation)
 1. [Extract and preprocess data from the MySQL database into a specified output directory](profiling)
 1. [Read preprocessed data from the specified output directory](reading)
 1. [Analysis and visualization](analysis)


# Analysis tools in R

Ensure that you always [update to the latest version](update)! Kindly
report any problems to [the admins](contact).

### Analysis 
[Clustering](clustering)  
[Core microbiota](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)  
[Correlation analysis](metrics)  
[Cross-hybridization](crosshyb)  
[Diversity analysis](diversity)  
[Functional Network Analysis](NetResponse)  
[Latent Class Analysis](LatentClassAnalysis)  
[Limma analysis](limma)  
[Group-wise comparisons](comparisons)  
[Phylogeny analysis](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)  
[RDA](rda)  
[ROC/AUC](roc)  
[Similarity/Dissimilarity matrices](metrics)  
[Stability analysis](stability)  

### Visualization
[Barplots](barplots)  
[Heatmaps](heatmap)  
[Visualization](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)  
[Density](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)  
[Hierarchical clustering](clustering)  
[Histograms, Clusters](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)  
[Interactive motionchart](https://github.com/microbiome/microbiome/blob/master/vignettes/Motionchart.Rmd)  

### Data preprocessing and I/O
[Extract data from SQL database (profiling)](https://github.com/microbiome/HITChipDB/blob/master/vignettes/vignette.Rmd)  
[Read preprocessed data](https://github.com/microbiome/microbiome/blob/master/vignettes/vignette.md)  
[SQL tools](https://github.com/microbiome/HITChipDB/blob/master/vignettes/vignette.Rmd)  

### Miscellaneous
[Experimental code](experimental)  


