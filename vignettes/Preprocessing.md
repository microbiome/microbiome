## Standard data processing operations


The high-quality [phyloseq package](http://joey711.github.io/phyloseq/) provides a complete set of tools for data subsetting, aggregating and filtering.

Download example data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling:


```r
library(microbiome)
pseq <- download_microbiome("dietswap")
```


### Sample operations

Sample names and variables


```r
head(sample_names(pseq))
```

```
## [1] "Sample-1" "Sample-2" "Sample-3" "Sample-4" "Sample-5" "Sample-6"
```

Sample sums


```r
head(sample_sums(pseq))
```

```
## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
##   533779  1330516  1822706   835998  1095023  1246234
```

Abundance of a given species in each sample


```r
head(get_sample(pseq, taxa_names(pseq)[1]))
```

```
## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
##       11       67       21       42       16       20
```

Filter samples


```r
f1 <- filterfun_sample(topp(0.1))
taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))
```

### Estimating relative abundancies

Estimate relative abundance of the taxa in each sample. Note: the
input data set needs to be in absolute scale (not logarithmic).


Calculating relative abundances for a phyloseq object:


```r
pseq <- transform_sample_counts(pseq, function(x) x/sum(x))
```


Calculating relative abundance for standard abundance matrix:


```r
dat <- otu_table(pseq)@.Data
rel <- relative.abundance(dat, det.th = 0)
```


Pick samples by specific metadata fields


```r
pseq.subset2 <- subset_samples(pseq.subset, nationality == "US")
```

```
## Error in sample_data(physeq): error in evaluating the argument 'object' in selecting a method for function 'sample_data': Error: object 'pseq.subset' not found
```


### Variable operations

Sample variable names


```r
sample_variables(pseq)
```

```
## [1] "subject"                "gender"                
## [3] "nationality"            "group"                 
## [5] "sample"                 "timepoint"             
## [7] "timepoint.within.group" "bmi_group"
```

Pick variable values for a given variable


```r
head(get_variable(pseq, sample_variables(pseq)[3]))
```

```
## [1] AAM AFR AFR AFR AFR AFR
## Levels: AAM AFR
```

```r
# Assign fields to sample metadata
# sample_data(GP)$human <- ..
```

### Taxa operations


Number of taxa


```r
ntaxa(pseq)
```

```
## [1] 130
```


Names


```r
rank_names(pseq)
```

```
## [1] "Phylum" "Genus"
```

```r
taxa_names(pseq)
```

```
##   [1] "Actinomycetaceae"                     
##   [2] "Aerococcus"                           
##   [3] "Aeromonas"                            
##   [4] "Akkermansia"                          
##   [5] "Alcaligenes faecalis et rel."         
##   [6] "Allistipes et rel."                   
##   [7] "Anaerobiospirillum"                   
##   [8] "Anaerofustis"                         
##   [9] "Anaerostipes caccae et rel."          
##  [10] "Anaerotruncus colihominis et rel."    
##  [11] "Anaerovorax odorimutans et rel."      
##  [12] "Aneurinibacillus"                     
##  [13] "Aquabacterium"                        
##  [14] "Asteroleplasma et rel."               
##  [15] "Atopobium"                            
##  [16] "Bacillus"                             
##  [17] "Bacteroides fragilis et rel."         
##  [18] "Bacteroides intestinalis et rel."     
##  [19] "Bacteroides ovatus et rel."           
##  [20] "Bacteroides plebeius et rel."         
##  [21] "Bacteroides splachnicus et rel."      
##  [22] "Bacteroides stercoris et rel."        
##  [23] "Bacteroides uniformis et rel."        
##  [24] "Bacteroides vulgatus et rel."         
##  [25] "Bifidobacterium"                      
##  [26] "Bilophila et rel."                    
##  [27] "Brachyspira"                          
##  [28] "Bryantella formatexigens et rel."     
##  [29] "Bulleidia moorei et rel."             
##  [30] "Burkholderia"                         
##  [31] "Butyrivibrio crossotus et rel."       
##  [32] "Campylobacter"                        
##  [33] "Catenibacterium mitsuokai et rel."    
##  [34] "Clostridium (sensu stricto)"          
##  [35] "Clostridium cellulosi et rel."        
##  [36] "Clostridium colinum et rel."          
##  [37] "Clostridium difficile et rel."        
##  [38] "Clostridium felsineum et rel."        
##  [39] "Clostridium leptum et rel."           
##  [40] "Clostridium nexile et rel."           
##  [41] "Clostridium orbiscindens et rel."     
##  [42] "Clostridium ramosum et rel."          
##  [43] "Clostridium sphenoides et rel."       
##  [44] "Clostridium stercorarium et rel."     
##  [45] "Clostridium symbiosum et rel."        
##  [46] "Clostridium thermocellum et rel."     
##  [47] "Collinsella"                          
##  [48] "Coprobacillus catenaformis et rel."   
##  [49] "Coprococcus eutactus et rel."         
##  [50] "Corynebacterium"                      
##  [51] "Desulfovibrio et rel."                
##  [52] "Dialister"                            
##  [53] "Dorea formicigenerans et rel."        
##  [54] "Eggerthella lenta et rel."            
##  [55] "Enterobacter aerogenes et rel."       
##  [56] "Enterococcus"                         
##  [57] "Escherichia coli et rel."             
##  [58] "Eubacterium biforme et rel."          
##  [59] "Eubacterium cylindroides et rel."     
##  [60] "Eubacterium hallii et rel."           
##  [61] "Eubacterium limosum et rel."          
##  [62] "Eubacterium rectale et rel."          
##  [63] "Eubacterium siraeum et rel."          
##  [64] "Eubacterium ventriosum et rel."       
##  [65] "Faecalibacterium prausnitzii et rel." 
##  [66] "Fusobacteria"                         
##  [67] "Gemella"                              
##  [68] "Granulicatella"                       
##  [69] "Haemophilus"                          
##  [70] "Helicobacter"                         
##  [71] "Klebisiella pneumoniae et rel."       
##  [72] "Lachnobacillus bovis et rel."         
##  [73] "Lachnospira pectinoschiza et rel."    
##  [74] "Lactobacillus catenaformis et rel."   
##  [75] "Lactobacillus gasseri et rel."        
##  [76] "Lactobacillus plantarum et rel."      
##  [77] "Lactobacillus salivarius et rel."     
##  [78] "Lactococcus"                          
##  [79] "Leminorella"                          
##  [80] "Megamonas hypermegale et rel."        
##  [81] "Megasphaera elsdenii et rel."         
##  [82] "Methylobacterium"                     
##  [83] "Micrococcaceae"                       
##  [84] "Mitsuokella multiacida et rel."       
##  [85] "Moraxellaceae"                        
##  [86] "Novosphingobium"                      
##  [87] "Oceanospirillum"                      
##  [88] "Oscillospira guillermondii et rel."   
##  [89] "Outgrouping clostridium cluster XIVa" 
##  [90] "Oxalobacter formigenes et rel."       
##  [91] "Papillibacter cinnamivorans et rel."  
##  [92] "Parabacteroides distasonis et rel."   
##  [93] "Peptococcus niger et rel."            
##  [94] "Peptostreptococcus anaerobius et rel."
##  [95] "Peptostreptococcus micros et rel."    
##  [96] "Phascolarctobacterium faecium et rel."
##  [97] "Prevotella melaninogenica et rel."    
##  [98] "Prevotella oralis et rel."            
##  [99] "Prevotella ruminicola et rel."        
## [100] "Prevotella tannerae et rel."          
## [101] "Propionibacterium"                    
## [102] "Proteus et rel."                      
## [103] "Pseudomonas"                          
## [104] "Roseburia intestinalis et rel."       
## [105] "Ruminococcus bromii et rel."          
## [106] "Ruminococcus callidus et rel."        
## [107] "Ruminococcus gnavus et rel."          
## [108] "Ruminococcus lactaris et rel."        
## [109] "Ruminococcus obeum et rel."           
## [110] "Serratia"                             
## [111] "Sporobacter termitidis et rel."       
## [112] "Staphylococcus"                       
## [113] "Streptococcus bovis et rel."          
## [114] "Streptococcus intermedius et rel."    
## [115] "Streptococcus mitis et rel."          
## [116] "Subdoligranulum variable at rel."     
## [117] "Sutterella wadsworthia et rel."       
## [118] "Tannerella et rel."                   
## [119] "Uncultured Bacteroidetes"             
## [120] "Uncultured Chroococcales"             
## [121] "Uncultured Clostridiales I"           
## [122] "Uncultured Clostridiales II"          
## [123] "Uncultured Mollicutes"                
## [124] "Uncultured Selenomonadaceae"          
## [125] "Veillonella"                          
## [126] "Vibrio"                               
## [127] "Weissella et rel."                    
## [128] "Wissella et rel."                     
## [129] "Xanthomonadaceae"                     
## [130] "Yersinia et rel."
```


Prune taxa


```r
taxa <- levelmap(NULL, "Phylum", "Genus", tax_table(pseq))$Bacteroidetes

# With given taxon names
ex2 <- prune_taxa(taxa, pseq)

# Taxa with positive sum across samples
ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)
```

Subset taxa


```r
pseq <- subset_taxa(pseq, Phylum == "Bacteroidetes")
```


Filter by user-specified function values (here variance)


```r
f <- filter_taxa(r, function(x) var(x) > 1e-05, TRUE)
```

```
## Error in inherits(x, get.component.classes()): object 'r' not found
```


List unique phylum-level groups: 


```r
head(get_taxa_unique(pseq, "Phylum"))
```

```
## [1] "Bacteroidetes"
```

Pick detected taxa by sample name:


```r
get_taxa(pseq, sample_names(pseq)[[1]])
```

```
##                 Allistipes et rel.       Bacteroides fragilis et rel. 
##                       0.0397580272                       0.0523156587 
##   Bacteroides intestinalis et rel.         Bacteroides ovatus et rel. 
##                       0.0014200634                       0.0504197430 
##       Bacteroides plebeius et rel.    Bacteroides splachnicus et rel. 
##                       0.0076304988                       0.0053111868 
##      Bacteroides stercoris et rel.      Bacteroides uniformis et rel. 
##                       0.0025235163                       0.0366855946 
##       Bacteroides vulgatus et rel. Parabacteroides distasonis et rel. 
##                       0.3279166097                       0.0220540711 
##  Prevotella melaninogenica et rel.          Prevotella oralis et rel. 
##                       0.0029281781                       0.0025965802 
##      Prevotella ruminicola et rel.        Prevotella tannerae et rel. 
##                       0.0001311404                       0.0074262944 
##                 Tannerella et rel.           Uncultured Bacteroidetes 
##                       0.0064596022                       0.0002116981
```


Taxa sums


```r
head(taxa_sums(pseq))
```

```
##               Allistipes et rel.     Bacteroides fragilis et rel. 
##                        4.8226896                        3.4909888 
## Bacteroides intestinalis et rel.       Bacteroides ovatus et rel. 
##                        0.3137315                        2.0882334 
##     Bacteroides plebeius et rel.  Bacteroides splachnicus et rel. 
##                        0.8411087                        1.0741668
```


### Merging operations for phyloseq objects


```r
merge_phyloseq()
merge_taxa()
tax_glom()
```


### Rarification


```r
pseq.rarified <- rarefy_even_depth(pseq)
```

```
## Error in validObject(.Object): invalid class "otu_table" object: 
##  OTU abundance data must have non-zero dimensions.
```


