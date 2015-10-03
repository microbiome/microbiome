### Prevalence of taxonomic groups


```r
# Load example data
library(microbiome)
pseq <- download_microbiome("peerj32")$physeq
```

```
## Downloading data set from Lahti et al. PeerJ, 2013: https://peerj.com/articles/32/
```

List prevalence measure for each group with a given detection threshold:


```r
head(prevalence(pseq, detection.threshold = 10, sort = FALSE))
```

```
##             Actinomycetaceae                   Aerococcus 
##                   0.13636364                   0.25000000 
##                    Aeromonas                  Akkermansia 
##                   0.31818182                   1.00000000 
## Alcaligenes faecalis et rel.           Allistipes et rel. 
##                   0.04545455                   1.00000000
```

List the taxa that are present over detection threshold in given
fraction of the samples:


```r
prevalent.taxa <- prevalent_taxa(pseq, detection.threshold = 50, prevalence.threshold = 0.2)
```


### Core microbiota

Determine core microbiota with the [blanket
analysis](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract)
based on various signal and prevalence thresholds.
 

```r
core <- core_matrix(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)))
```

### Core 2D line plots


```r
p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 10^(0:4)), plot.type = "lineplot")
```

![plot of chunk core-example2](figure/core-example2-1.png) 

### Core heatmaps


```r
p <- plot_core(pseq, prevalence.intervals = seq(10, 100, 10), detection.thresholds = c(0, 2^(0:14)), plot.type = "heatmap")
```

![plot of chunk core-example3](figure/core-example3-1.png) 



```r
p <- core_heatmap(pseq, palette = "spectral")
print(p)
```

```
## $plot
```

![plot of chunk core-example4](figure/core-example4-1.png) 

```
## 
## $data
##                                            0.001 0.00134581030387798
## Actinomycetaceae                        0.000000            0.000000
## Aerococcus                              0.000000            0.000000
## Aeromonas                               4.545455            2.272727
## Akkermansia                           100.000000          100.000000
## Alcaligenes faecalis et rel.            0.000000            0.000000
## Allistipes et rel.                     97.727273           97.727273
## Anaerobiospirillum                      0.000000            0.000000
## Anaerofustis                           45.454545           25.000000
## Anaerostipes caccae et rel.           100.000000          100.000000
## Anaerotruncus colihominis et rel.      84.090909           77.272727
## Anaerovorax odorimutans et rel.        97.727273           93.181818
## Aneurinibacillus                       22.727273           13.636364
## Aquabacterium                          56.818182           50.000000
## Asteroleplasma et rel.                  0.000000            0.000000
## Atopobium                              20.454545            6.818182
## Bacillus                                9.090909            0.000000
## Bacteroides fragilis et rel.          100.000000           97.727273
## Bacteroides intestinalis et rel.      100.000000          100.000000
## Bacteroides ovatus et rel.            100.000000           97.727273
## Bacteroides plebeius et rel.          100.000000           97.727273
## Bacteroides splachnicus et rel.        97.727273           93.181818
## Bacteroides stercoris et rel.          90.909091           79.545455
## Bacteroides uniformis et rel.         100.000000          100.000000
## Bacteroides vulgatus et rel.          100.000000          100.000000
## Bifidobacterium                       100.000000           97.727273
## Bilophila et rel.                       2.272727            2.272727
## Brachyspira                             0.000000            0.000000
## Bryantella formatexigens et rel.      100.000000          100.000000
## Bulleidia moorei et rel.               13.636364           11.363636
## Burkholderia                           52.272727           36.363636
## Butyrivibrio crossotus et rel.        100.000000          100.000000
## Campylobacter                           0.000000            0.000000
## Catenibacterium mitsuokai et rel.      18.181818           18.181818
## Clostridium (sensu stricto)             9.090909            6.818182
## Clostridium cellulosi et rel.         100.000000           97.727273
## Clostridium colinum et rel.            95.454545           86.363636
## Clostridium difficile et rel.          77.272727           59.090909
## Clostridium felsineum et rel.           4.545455            0.000000
## Clostridium leptum et rel.            100.000000          100.000000
## Clostridium nexile et rel.            100.000000          100.000000
## Clostridium orbiscindens et rel.      100.000000          100.000000
## Clostridium ramosum et rel.            79.545455           72.727273
## Clostridium sphenoides et rel.        100.000000          100.000000
## Clostridium stercorarium et rel.       84.090909           77.272727
## Clostridium symbiosum et rel.         100.000000          100.000000
## Clostridium thermocellum et rel.        0.000000            0.000000
## Collinsella                            97.727273           84.090909
## Coprobacillus catenaformis et rel.     77.272727           75.000000
## Coprococcus eutactus et rel.          100.000000          100.000000
## Corynebacterium                         0.000000            0.000000
## Desulfovibrio et rel.                  13.636364           11.363636
## Dialister                              40.909091           38.636364
## Dorea formicigenerans et rel.         100.000000          100.000000
## Eggerthella lenta et rel.              59.090909           38.636364
## Enterobacter aerogenes et rel.         11.363636            6.818182
## Enterococcus                           13.636364           11.363636
## Escherichia coli et rel.                4.545455            4.545455
## Eubacterium biforme et rel.            77.272727           65.909091
## Eubacterium cylindroides et rel.       18.181818           15.909091
## Eubacterium hallii et rel.            100.000000          100.000000
## Eubacterium limosum et rel.             6.818182            0.000000
## Eubacterium rectale et rel.           100.000000          100.000000
## Eubacterium siraeum et rel.            88.636364           77.272727
## Eubacterium ventriosum et rel.        100.000000          100.000000
## Faecalibacterium prausnitzii et rel.  100.000000          100.000000
## Fusobacteria                            0.000000            0.000000
## Gemella                                 0.000000            0.000000
## Granulicatella                         20.454545           13.636364
## Haemophilus                             9.090909            6.818182
## Helicobacter                            0.000000            0.000000
## Klebisiella pneumoniae et rel.          9.090909            2.272727
## Lachnobacillus bovis et rel.          100.000000          100.000000
## Lachnospira pectinoschiza et rel.     100.000000          100.000000
## Lactobacillus catenaformis et rel.      0.000000            0.000000
## Lactobacillus gasseri et rel.           4.545455            4.545455
## Lactobacillus plantarum et rel.         6.818182            2.272727
## Lactobacillus salivarius et rel.        2.272727            0.000000
## Lactococcus                            40.909091           31.818182
## Leminorella                             4.545455            2.272727
## Megamonas hypermegale et rel.           4.545455            2.272727
## Megasphaera elsdenii et rel.            9.090909            9.090909
## Methylobacterium                        0.000000            0.000000
## Micrococcaceae                          0.000000            0.000000
## Mitsuokella multiacida et rel.          4.545455            4.545455
## Moraxellaceae                           0.000000            0.000000
## Novosphingobium                        18.181818           13.636364
## Oceanospirillum                         4.545455            2.272727
## Oscillospira guillermondii et rel.    100.000000           97.727273
## Outgrouping clostridium cluster XIVa  100.000000          100.000000
## Oxalobacter formigenes et rel.        100.000000           97.727273
## Papillibacter cinnamivorans et rel.   100.000000          100.000000
## Parabacteroides distasonis et rel.    100.000000          100.000000
## Peptococcus niger et rel.              45.454545           34.090909
## Peptostreptococcus anaerobius et rel.   0.000000            0.000000
## Peptostreptococcus micros et rel.       0.000000            0.000000
## Phascolarctobacterium faecium et rel.   9.090909            4.545455
## Prevotella melaninogenica et rel.      29.545455           22.727273
## Prevotella oralis et rel.              56.818182           43.181818
## Prevotella ruminicola et rel.          86.363636           77.272727
## Prevotella tannerae et rel.           100.000000          100.000000
## Propionibacterium                      13.636364            6.818182
## Proteus et rel.                         0.000000            0.000000
## Pseudomonas                            13.636364            9.090909
## Roseburia intestinalis et rel.        100.000000          100.000000
## Ruminococcus bromii et rel.            95.454545           95.454545
## Ruminococcus callidus et rel.         100.000000          100.000000
## Ruminococcus gnavus et rel.           100.000000          100.000000
## Ruminococcus lactaris et rel.         100.000000          100.000000
## Ruminococcus obeum et rel.            100.000000          100.000000
## Serratia                               13.636364            9.090909
## Sporobacter termitidis et rel.        100.000000          100.000000
## Staphylococcus                          4.545455            4.545455
## Streptococcus bovis et rel.           100.000000          100.000000
## Streptococcus intermedius et rel.     100.000000           97.727273
## Streptococcus mitis et rel.           100.000000          100.000000
## Subdoligranulum variable at rel.      100.000000          100.000000
## Sutterella wadsworthia et rel.         61.363636           47.727273
## Tannerella et rel.                     93.181818           88.636364
## Uncultured Bacteroidetes               77.272727           65.909091
## Uncultured Chroococcales               18.181818           15.909091
## Uncultured Clostridiales I             56.818182           54.545455
## Uncultured Clostridiales II            47.727273           36.363636
## Uncultured Mollicutes                  34.090909           29.545455
## Uncultured Selenomonadaceae            13.636364           11.363636
## Veillonella                            34.090909           27.272727
## Weissella et rel.                      15.909091           13.636364
## Vibrio                                  6.818182            2.272727
## Wissella et rel.                       11.363636            6.818182
## Xanthomonadaceae                        9.090909            4.545455
## Yersinia et rel.                        0.000000            0.000000
##                                       0.00181120537402416
## Actinomycetaceae                                 0.000000
## Aerococcus                                       0.000000
## Aeromonas                                        0.000000
## Akkermansia                                    100.000000
## Alcaligenes faecalis et rel.                     0.000000
## Allistipes et rel.                              93.181818
## Anaerobiospirillum                               0.000000
## Anaerofustis                                    13.636364
## Anaerostipes caccae et rel.                    100.000000
## Anaerotruncus colihominis et rel.               72.727273
## Anaerovorax odorimutans et rel.                 90.909091
## Aneurinibacillus                                 6.818182
## Aquabacterium                                   36.363636
## Asteroleplasma et rel.                           0.000000
## Atopobium                                        4.545455
## Bacillus                                         0.000000
## Bacteroides fragilis et rel.                    97.727273
## Bacteroides intestinalis et rel.                97.727273
## Bacteroides ovatus et rel.                      88.636364
## Bacteroides plebeius et rel.                    97.727273
## Bacteroides splachnicus et rel.                 90.909091
## Bacteroides stercoris et rel.                   72.727273
## Bacteroides uniformis et rel.                  100.000000
## Bacteroides vulgatus et rel.                   100.000000
## Bifidobacterium                                 97.727273
## Bilophila et rel.                                2.272727
## Brachyspira                                      0.000000
## Bryantella formatexigens et rel.               100.000000
## Bulleidia moorei et rel.                         6.818182
## Burkholderia                                    25.000000
## Butyrivibrio crossotus et rel.                 100.000000
## Campylobacter                                    0.000000
## Catenibacterium mitsuokai et rel.               13.636364
## Clostridium (sensu stricto)                      2.272727
## Clostridium cellulosi et rel.                   97.727273
## Clostridium colinum et rel.                     70.454545
## Clostridium difficile et rel.                   38.636364
## Clostridium felsineum et rel.                    0.000000
## Clostridium leptum et rel.                     100.000000
## Clostridium nexile et rel.                     100.000000
## Clostridium orbiscindens et rel.               100.000000
## Clostridium ramosum et rel.                     54.545455
## Clostridium sphenoides et rel.                 100.000000
## Clostridium stercorarium et rel.                59.090909
## Clostridium symbiosum et rel.                  100.000000
## Clostridium thermocellum et rel.                 0.000000
## Collinsella                                     75.000000
## Coprobacillus catenaformis et rel.              59.090909
## Coprococcus eutactus et rel.                   100.000000
## Corynebacterium                                  0.000000
## Desulfovibrio et rel.                            2.272727
## Dialister                                       38.636364
## Dorea formicigenerans et rel.                  100.000000
## Eggerthella lenta et rel.                       15.909091
## Enterobacter aerogenes et rel.                   4.545455
## Enterococcus                                     4.545455
## Escherichia coli et rel.                         4.545455
## Eubacterium biforme et rel.                     52.272727
## Eubacterium cylindroides et rel.                13.636364
## Eubacterium hallii et rel.                     100.000000
## Eubacterium limosum et rel.                      0.000000
## Eubacterium rectale et rel.                    100.000000
## Eubacterium siraeum et rel.                     50.000000
## Eubacterium ventriosum et rel.                 100.000000
## Faecalibacterium prausnitzii et rel.           100.000000
## Fusobacteria                                     0.000000
## Gemella                                          0.000000
## Granulicatella                                   9.090909
## Haemophilus                                      4.545455
## Helicobacter                                     0.000000
## Klebisiella pneumoniae et rel.                   2.272727
## Lachnobacillus bovis et rel.                   100.000000
## Lachnospira pectinoschiza et rel.              100.000000
## Lactobacillus catenaformis et rel.               0.000000
## Lactobacillus gasseri et rel.                    4.545455
## Lactobacillus plantarum et rel.                  0.000000
## Lactobacillus salivarius et rel.                 0.000000
## Lactococcus                                     25.000000
## Leminorella                                      0.000000
## Megamonas hypermegale et rel.                    2.272727
## Megasphaera elsdenii et rel.                     4.545455
## Methylobacterium                                 0.000000
## Micrococcaceae                                   0.000000
## Mitsuokella multiacida et rel.                   4.545455
## Moraxellaceae                                    0.000000
## Novosphingobium                                  6.818182
## Oceanospirillum                                  0.000000
## Oscillospira guillermondii et rel.              97.727273
## Outgrouping clostridium cluster XIVa           100.000000
## Oxalobacter formigenes et rel.                  95.454545
## Papillibacter cinnamivorans et rel.            100.000000
## Parabacteroides distasonis et rel.              93.181818
## Peptococcus niger et rel.                       22.727273
## Peptostreptococcus anaerobius et rel.            0.000000
## Peptostreptococcus micros et rel.                0.000000
## Phascolarctobacterium faecium et rel.            2.272727
## Prevotella melaninogenica et rel.               22.727273
## Prevotella oralis et rel.                       31.818182
## Prevotella ruminicola et rel.                   50.000000
## Prevotella tannerae et rel.                     95.454545
## Propionibacterium                                4.545455
## Proteus et rel.                                  0.000000
## Pseudomonas                                      4.545455
## Roseburia intestinalis et rel.                 100.000000
## Ruminococcus bromii et rel.                     88.636364
## Ruminococcus callidus et rel.                   88.636364
## Ruminococcus gnavus et rel.                    100.000000
## Ruminococcus lactaris et rel.                  100.000000
## Ruminococcus obeum et rel.                     100.000000
## Serratia                                         6.818182
## Sporobacter termitidis et rel.                  97.727273
## Staphylococcus                                   0.000000
## Streptococcus bovis et rel.                    100.000000
## Streptococcus intermedius et rel.               95.454545
## Streptococcus mitis et rel.                    100.000000
## Subdoligranulum variable at rel.               100.000000
## Sutterella wadsworthia et rel.                  38.636364
## Tannerella et rel.                              75.000000
## Uncultured Bacteroidetes                        56.818182
## Uncultured Chroococcales                         9.090909
## Uncultured Clostridiales I                      47.727273
## Uncultured Clostridiales II                     25.000000
## Uncultured Mollicutes                           18.181818
## Uncultured Selenomonadaceae                      2.272727
## Veillonella                                     22.727273
## Weissella et rel.                                9.090909
## Vibrio                                           2.272727
## Wissella et rel.                                 2.272727
## Xanthomonadaceae                                 2.272727
## Yersinia et rel.                                 0.000000
##                                       0.00243753885480089
## Actinomycetaceae                                 0.000000
## Aerococcus                                       0.000000
## Aeromonas                                        0.000000
## Akkermansia                                    100.000000
## Alcaligenes faecalis et rel.                     0.000000
## Allistipes et rel.                              90.909091
## Anaerobiospirillum                               0.000000
## Anaerofustis                                     9.090909
## Anaerostipes caccae et rel.                    100.000000
## Anaerotruncus colihominis et rel.               50.000000
## Anaerovorax odorimutans et rel.                 79.545455
## Aneurinibacillus                                 2.272727
## Aquabacterium                                   27.272727
## Asteroleplasma et rel.                           0.000000
## Atopobium                                        2.272727
## Bacillus                                         0.000000
## Bacteroides fragilis et rel.                    90.909091
## Bacteroides intestinalis et rel.                97.727273
## Bacteroides ovatus et rel.                      77.272727
## Bacteroides plebeius et rel.                    93.181818
## Bacteroides splachnicus et rel.                 79.545455
## Bacteroides stercoris et rel.                   56.818182
## Bacteroides uniformis et rel.                  100.000000
## Bacteroides vulgatus et rel.                   100.000000
## Bifidobacterium                                 93.181818
## Bilophila et rel.                                0.000000
## Brachyspira                                      0.000000
## Bryantella formatexigens et rel.                97.727273
## Bulleidia moorei et rel.                         0.000000
## Burkholderia                                     9.090909
## Butyrivibrio crossotus et rel.                 100.000000
## Campylobacter                                    0.000000
## Catenibacterium mitsuokai et rel.                9.090909
## Clostridium (sensu stricto)                      0.000000
## Clostridium cellulosi et rel.                   93.181818
## Clostridium colinum et rel.                     59.090909
## Clostridium difficile et rel.                   27.272727
## Clostridium felsineum et rel.                    0.000000
## Clostridium leptum et rel.                      97.727273
## Clostridium nexile et rel.                     100.000000
## Clostridium orbiscindens et rel.                97.727273
## Clostridium ramosum et rel.                     36.363636
## Clostridium sphenoides et rel.                 100.000000
## Clostridium stercorarium et rel.                47.727273
## Clostridium symbiosum et rel.                  100.000000
## Clostridium thermocellum et rel.                 0.000000
## Collinsella                                     59.090909
## Coprobacillus catenaformis et rel.              40.909091
## Coprococcus eutactus et rel.                   100.000000
## Corynebacterium                                  0.000000
## Desulfovibrio et rel.                            0.000000
## Dialister                                       34.090909
## Dorea formicigenerans et rel.                  100.000000
## Eggerthella lenta et rel.                        6.818182
## Enterobacter aerogenes et rel.                   2.272727
## Enterococcus                                     0.000000
## Escherichia coli et rel.                         2.272727
## Eubacterium biforme et rel.                     34.090909
## Eubacterium cylindroides et rel.                 9.090909
## Eubacterium hallii et rel.                     100.000000
## Eubacterium limosum et rel.                      0.000000
## Eubacterium rectale et rel.                    100.000000
## Eubacterium siraeum et rel.                     25.000000
## Eubacterium ventriosum et rel.                 100.000000
## Faecalibacterium prausnitzii et rel.           100.000000
## Fusobacteria                                     0.000000
## Gemella                                          0.000000
## Granulicatella                                   2.272727
## Haemophilus                                      4.545455
## Helicobacter                                     0.000000
## Klebisiella pneumoniae et rel.                   2.272727
## Lachnobacillus bovis et rel.                   100.000000
## Lachnospira pectinoschiza et rel.              100.000000
## Lactobacillus catenaformis et rel.               0.000000
## Lactobacillus gasseri et rel.                    2.272727
## Lactobacillus plantarum et rel.                  0.000000
## Lactobacillus salivarius et rel.                 0.000000
## Lactococcus                                     25.000000
## Leminorella                                      0.000000
## Megamonas hypermegale et rel.                    0.000000
## Megasphaera elsdenii et rel.                     4.545455
## Methylobacterium                                 0.000000
## Micrococcaceae                                   0.000000
## Mitsuokella multiacida et rel.                   4.545455
## Moraxellaceae                                    0.000000
## Novosphingobium                                  2.272727
## Oceanospirillum                                  0.000000
## Oscillospira guillermondii et rel.              97.727273
## Outgrouping clostridium cluster XIVa           100.000000
## Oxalobacter formigenes et rel.                  93.181818
## Papillibacter cinnamivorans et rel.            100.000000
## Parabacteroides distasonis et rel.              90.909091
## Peptococcus niger et rel.                        6.818182
## Peptostreptococcus anaerobius et rel.            0.000000
## Peptostreptococcus micros et rel.                0.000000
## Phascolarctobacterium faecium et rel.            0.000000
## Prevotella melaninogenica et rel.               22.727273
## Prevotella oralis et rel.                       18.181818
## Prevotella ruminicola et rel.                   27.272727
## Prevotella tannerae et rel.                     93.181818
## Propionibacterium                                2.272727
## Proteus et rel.                                  0.000000
## Pseudomonas                                      4.545455
## Roseburia intestinalis et rel.                 100.000000
## Ruminococcus bromii et rel.                     86.363636
## Ruminococcus callidus et rel.                   77.272727
## Ruminococcus gnavus et rel.                    100.000000
## Ruminococcus lactaris et rel.                  100.000000
## Ruminococcus obeum et rel.                     100.000000
## Serratia                                         6.818182
## Sporobacter termitidis et rel.                  88.636364
## Staphylococcus                                   0.000000
## Streptococcus bovis et rel.                    100.000000
## Streptococcus intermedius et rel.               79.545455
## Streptococcus mitis et rel.                    100.000000
## Subdoligranulum variable at rel.               100.000000
## Sutterella wadsworthia et rel.                  22.727273
## Tannerella et rel.                              54.545455
## Uncultured Bacteroidetes                        36.363636
## Uncultured Chroococcales                         9.090909
## Uncultured Clostridiales I                      36.363636
## Uncultured Clostridiales II                     11.363636
## Uncultured Mollicutes                           15.909091
## Uncultured Selenomonadaceae                      2.272727
## Veillonella                                     18.181818
## Weissella et rel.                                2.272727
## Vibrio                                           0.000000
## Wissella et rel.                                 0.000000
## Xanthomonadaceae                                 0.000000
## Yersinia et rel.                                 0.000000
##                                       0.00328046490689398
## Actinomycetaceae                                 0.000000
## Aerococcus                                       0.000000
## Aeromonas                                        0.000000
## Akkermansia                                    100.000000
## Alcaligenes faecalis et rel.                     0.000000
## Allistipes et rel.                              88.636364
## Anaerobiospirillum                               0.000000
## Anaerofustis                                     6.818182
## Anaerostipes caccae et rel.                    100.000000
## Anaerotruncus colihominis et rel.               29.545455
## Anaerovorax odorimutans et rel.                 59.090909
## Aneurinibacillus                                 0.000000
## Aquabacterium                                   11.363636
## Asteroleplasma et rel.                           0.000000
## Atopobium                                        0.000000
## Bacillus                                         0.000000
## Bacteroides fragilis et rel.                    84.090909
## Bacteroides intestinalis et rel.                97.727273
## Bacteroides ovatus et rel.                      59.090909
## Bacteroides plebeius et rel.                    88.636364
## Bacteroides splachnicus et rel.                 50.000000
## Bacteroides stercoris et rel.                   50.000000
## Bacteroides uniformis et rel.                   97.727273
## Bacteroides vulgatus et rel.                   100.000000
## Bifidobacterium                                 90.909091
## Bilophila et rel.                                0.000000
## Brachyspira                                      0.000000
## Bryantella formatexigens et rel.                95.454545
## Bulleidia moorei et rel.                         0.000000
## Burkholderia                                     2.272727
## Butyrivibrio crossotus et rel.                 100.000000
## Campylobacter                                    0.000000
## Catenibacterium mitsuokai et rel.                9.090909
## Clostridium (sensu stricto)                      0.000000
## Clostridium cellulosi et rel.                   86.363636
## Clostridium colinum et rel.                     47.727273
## Clostridium difficile et rel.                    9.090909
## Clostridium felsineum et rel.                    0.000000
## Clostridium leptum et rel.                      93.181818
## Clostridium nexile et rel.                     100.000000
## Clostridium orbiscindens et rel.                95.454545
## Clostridium ramosum et rel.                     18.181818
## Clostridium sphenoides et rel.                  97.727273
## Clostridium stercorarium et rel.                13.636364
## Clostridium symbiosum et rel.                  100.000000
## Clostridium thermocellum et rel.                 0.000000
## Collinsella                                     52.272727
## Coprobacillus catenaformis et rel.              29.545455
## Coprococcus eutactus et rel.                   100.000000
## Corynebacterium                                  0.000000
## Desulfovibrio et rel.                            0.000000
## Dialister                                       27.272727
## Dorea formicigenerans et rel.                  100.000000
## Eggerthella lenta et rel.                        0.000000
## Enterobacter aerogenes et rel.                   0.000000
## Enterococcus                                     0.000000
## Escherichia coli et rel.                         2.272727
## Eubacterium biforme et rel.                     22.727273
## Eubacterium cylindroides et rel.                 4.545455
## Eubacterium hallii et rel.                     100.000000
## Eubacterium limosum et rel.                      0.000000
## Eubacterium rectale et rel.                    100.000000
## Eubacterium siraeum et rel.                     15.909091
## Eubacterium ventriosum et rel.                 100.000000
## Faecalibacterium prausnitzii et rel.           100.000000
## Fusobacteria                                     0.000000
## Gemella                                          0.000000
## Granulicatella                                   0.000000
## Haemophilus                                      4.545455
## Helicobacter                                     0.000000
## Klebisiella pneumoniae et rel.                   2.272727
## Lachnobacillus bovis et rel.                   100.000000
## Lachnospira pectinoschiza et rel.              100.000000
## Lactobacillus catenaformis et rel.               0.000000
## Lactobacillus gasseri et rel.                    2.272727
## Lactobacillus plantarum et rel.                  0.000000
## Lactobacillus salivarius et rel.                 0.000000
## Lactococcus                                     20.454545
## Leminorella                                      0.000000
## Megamonas hypermegale et rel.                    0.000000
## Megasphaera elsdenii et rel.                     4.545455
## Methylobacterium                                 0.000000
## Micrococcaceae                                   0.000000
## Mitsuokella multiacida et rel.                   4.545455
## Moraxellaceae                                    0.000000
## Novosphingobium                                  0.000000
## Oceanospirillum                                  0.000000
## Oscillospira guillermondii et rel.              97.727273
## Outgrouping clostridium cluster XIVa           100.000000
## Oxalobacter formigenes et rel.                  75.000000
## Papillibacter cinnamivorans et rel.            100.000000
## Parabacteroides distasonis et rel.              81.818182
## Peptococcus niger et rel.                        2.272727
## Peptostreptococcus anaerobius et rel.            0.000000
## Peptostreptococcus micros et rel.                0.000000
## Phascolarctobacterium faecium et rel.            0.000000
## Prevotella melaninogenica et rel.               18.181818
## Prevotella oralis et rel.                       15.909091
## Prevotella ruminicola et rel.                   11.363636
## Prevotella tannerae et rel.                     81.818182
## Propionibacterium                                0.000000
## Proteus et rel.                                  0.000000
## Pseudomonas                                      4.545455
## Roseburia intestinalis et rel.                 100.000000
## Ruminococcus bromii et rel.                     86.363636
## Ruminococcus callidus et rel.                   68.181818
## Ruminococcus gnavus et rel.                    100.000000
## Ruminococcus lactaris et rel.                  100.000000
## Ruminococcus obeum et rel.                     100.000000
## Serratia                                         6.818182
## Sporobacter termitidis et rel.                  75.000000
## Staphylococcus                                   0.000000
## Streptococcus bovis et rel.                    100.000000
## Streptococcus intermedius et rel.               70.454545
## Streptococcus mitis et rel.                    100.000000
## Subdoligranulum variable at rel.               100.000000
## Sutterella wadsworthia et rel.                   6.818182
## Tannerella et rel.                              31.818182
## Uncultured Bacteroidetes                        31.818182
## Uncultured Chroococcales                         6.818182
## Uncultured Clostridiales I                      31.818182
## Uncultured Clostridiales II                      9.090909
## Uncultured Mollicutes                            9.090909
## Uncultured Selenomonadaceae                      2.272727
## Veillonella                                     11.363636
## Weissella et rel.                                0.000000
## Vibrio                                           0.000000
## Wissella et rel.                                 0.000000
## Xanthomonadaceae                                 0.000000
## Yersinia et rel.                                 0.000000
##                                       0.00441488347320805
## Actinomycetaceae                                 0.000000
## Aerococcus                                       0.000000
## Aeromonas                                        0.000000
## Akkermansia                                     95.454545
## Alcaligenes faecalis et rel.                     0.000000
## Allistipes et rel.                              75.000000
## Anaerobiospirillum                               0.000000
## Anaerofustis                                     4.545455
## Anaerostipes caccae et rel.                    100.000000
## Anaerotruncus colihominis et rel.               18.181818
## Anaerovorax odorimutans et rel.                 29.545455
## Aneurinibacillus                                 0.000000
## Aquabacterium                                    2.272727
## Asteroleplasma et rel.                           0.000000
## Atopobium                                        0.000000
## Bacillus                                         0.000000
## Bacteroides fragilis et rel.                    68.181818
## Bacteroides intestinalis et rel.                86.363636
## Bacteroides ovatus et rel.                      47.727273
## Bacteroides plebeius et rel.                    68.181818
## Bacteroides splachnicus et rel.                 27.272727
## Bacteroides stercoris et rel.                   36.363636
## Bacteroides uniformis et rel.                   88.636364
## Bacteroides vulgatus et rel.                    97.727273
## Bifidobacterium                                 90.909091
## Bilophila et rel.                                0.000000
## Brachyspira                                      0.000000
## Bryantella formatexigens et rel.                95.454545
## Bulleidia moorei et rel.                         0.000000
## Burkholderia                                     0.000000
## Butyrivibrio crossotus et rel.                  88.636364
## Campylobacter                                    0.000000
## Catenibacterium mitsuokai et rel.                6.818182
## Clostridium (sensu stricto)                      0.000000
## Clostridium cellulosi et rel.                   65.909091
## Clostridium colinum et rel.                     15.909091
## Clostridium difficile et rel.                    4.545455
## Clostridium felsineum et rel.                    0.000000
## Clostridium leptum et rel.                      90.909091
## Clostridium nexile et rel.                     100.000000
## Clostridium orbiscindens et rel.                93.181818
## Clostridium ramosum et rel.                      6.818182
## Clostridium sphenoides et rel.                  95.454545
## Clostridium stercorarium et rel.                 4.545455
## Clostridium symbiosum et rel.                  100.000000
## Clostridium thermocellum et rel.                 0.000000
## Collinsella                                     38.636364
## Coprobacillus catenaformis et rel.              13.636364
## Coprococcus eutactus et rel.                   100.000000
## Corynebacterium                                  0.000000
## Desulfovibrio et rel.                            0.000000
## Dialister                                       20.454545
## Dorea formicigenerans et rel.                  100.000000
## Eggerthella lenta et rel.                        0.000000
## Enterobacter aerogenes et rel.                   0.000000
## Enterococcus                                     0.000000
## Escherichia coli et rel.                         2.272727
## Eubacterium biforme et rel.                     18.181818
## Eubacterium cylindroides et rel.                 4.545455
## Eubacterium hallii et rel.                     100.000000
## Eubacterium limosum et rel.                      0.000000
## Eubacterium rectale et rel.                    100.000000
## Eubacterium siraeum et rel.                      9.090909
## Eubacterium ventriosum et rel.                 100.000000
## Faecalibacterium prausnitzii et rel.           100.000000
## Fusobacteria                                     0.000000
## Gemella                                          0.000000
## Granulicatella                                   0.000000
## Haemophilus                                      4.545455
## Helicobacter                                     0.000000
## Klebisiella pneumoniae et rel.                   2.272727
## Lachnobacillus bovis et rel.                   100.000000
## Lachnospira pectinoschiza et rel.              100.000000
## Lactobacillus catenaformis et rel.               0.000000
## Lactobacillus gasseri et rel.                    2.272727
## Lactobacillus plantarum et rel.                  0.000000
## Lactobacillus salivarius et rel.                 0.000000
## Lactococcus                                     20.454545
## Leminorella                                      0.000000
## Megamonas hypermegale et rel.                    0.000000
## Megasphaera elsdenii et rel.                     0.000000
## Methylobacterium                                 0.000000
## Micrococcaceae                                   0.000000
## Mitsuokella multiacida et rel.                   4.545455
## Moraxellaceae                                    0.000000
## Novosphingobium                                  0.000000
## Oceanospirillum                                  0.000000
## Oscillospira guillermondii et rel.              93.181818
## Outgrouping clostridium cluster XIVa            95.454545
## Oxalobacter formigenes et rel.                  63.636364
## Papillibacter cinnamivorans et rel.            100.000000
## Parabacteroides distasonis et rel.              63.636364
## Peptococcus niger et rel.                        0.000000
## Peptostreptococcus anaerobius et rel.            0.000000
## Peptostreptococcus micros et rel.                0.000000
## Phascolarctobacterium faecium et rel.            0.000000
## Prevotella melaninogenica et rel.               13.636364
## Prevotella oralis et rel.                       11.363636
## Prevotella ruminicola et rel.                    6.818182
## Prevotella tannerae et rel.                     61.363636
## Propionibacterium                                0.000000
## Proteus et rel.                                  0.000000
## Pseudomonas                                      2.272727
## Roseburia intestinalis et rel.                 100.000000
## Ruminococcus bromii et rel.                     86.363636
## Ruminococcus callidus et rel.                   52.272727
## Ruminococcus gnavus et rel.                    100.000000
## Ruminococcus lactaris et rel.                  100.000000
## Ruminococcus obeum et rel.                     100.000000
## Serratia                                         4.545455
## Sporobacter termitidis et rel.                  54.545455
## Staphylococcus                                   0.000000
## Streptococcus bovis et rel.                     97.727273
## Streptococcus intermedius et rel.               45.454545
## Streptococcus mitis et rel.                     95.454545
## Subdoligranulum variable at rel.                97.727273
## Sutterella wadsworthia et rel.                   0.000000
## Tannerella et rel.                              18.181818
## Uncultured Bacteroidetes                        20.454545
## Uncultured Chroococcales                         2.272727
## Uncultured Clostridiales I                      27.272727
## Uncultured Clostridiales II                      0.000000
## Uncultured Mollicutes                            2.272727
## Uncultured Selenomonadaceae                      2.272727
## Veillonella                                     11.363636
## Weissella et rel.                                0.000000
## Vibrio                                           0.000000
## Wissella et rel.                                 0.000000
## Xanthomonadaceae                                 0.000000
## Yersinia et rel.                                 0.000000
##                                       0.00594159566866402
## Actinomycetaceae                                 0.000000
## Aerococcus                                       0.000000
## Aeromonas                                        0.000000
## Akkermansia                                     84.090909
## Alcaligenes faecalis et rel.                     0.000000
## Allistipes et rel.                              68.181818
## Anaerobiospirillum                               0.000000
## Anaerofustis                                     4.545455
## Anaerostipes caccae et rel.                    100.000000
## Anaerotruncus colihominis et rel.                4.545455
## Anaerovorax odorimutans et rel.                 18.181818
## Aneurinibacillus                                 0.000000
## Aquabacterium                                    2.272727
## Asteroleplasma et rel.                           0.000000
## Atopobium                                        0.000000
## Bacillus                                         0.000000
## Bacteroides fragilis et rel.                    47.727273
## Bacteroides intestinalis et rel.                79.545455
## Bacteroides ovatus et rel.                      31.818182
## Bacteroides plebeius et rel.                    36.363636
## Bacteroides splachnicus et rel.                  9.090909
## Bacteroides stercoris et rel.                   29.545455
## Bacteroides uniformis et rel.                   81.818182
## Bacteroides vulgatus et rel.                    95.454545
## Bifidobacterium                                 88.636364
## Bilophila et rel.                                0.000000
## Brachyspira                                      0.000000
## Bryantella formatexigens et rel.                95.454545
## Bulleidia moorei et rel.                         0.000000
## Burkholderia                                     0.000000
## Butyrivibrio crossotus et rel.                  52.272727
## Campylobacter                                    0.000000
## Catenibacterium mitsuokai et rel.                6.818182
## Clostridium (sensu stricto)                      0.000000
## Clostridium cellulosi et rel.                   54.545455
## Clostridium colinum et rel.                     11.363636
## Clostridium difficile et rel.                    2.272727
## Clostridium felsineum et rel.                    0.000000
## Clostridium leptum et rel.                      72.727273
## Clostridium nexile et rel.                     100.000000
## Clostridium orbiscindens et rel.                75.000000
## Clostridium ramosum et rel.                      2.272727
## Clostridium sphenoides et rel.                  95.454545
## Clostridium stercorarium et rel.                 0.000000
## Clostridium symbiosum et rel.                   97.727273
## Clostridium thermocellum et rel.                 0.000000
## Collinsella                                     20.454545
## Coprobacillus catenaformis et rel.               4.545455
## Coprococcus eutactus et rel.                   100.000000
## Corynebacterium                                  0.000000
## Desulfovibrio et rel.                            0.000000
## Dialister                                       15.909091
## Dorea formicigenerans et rel.                  100.000000
## Eggerthella lenta et rel.                        0.000000
## Enterobacter aerogenes et rel.                   0.000000
## Enterococcus                                     0.000000
## Escherichia coli et rel.                         2.272727
## Eubacterium biforme et rel.                     15.909091
## Eubacterium cylindroides et rel.                 4.545455
## Eubacterium hallii et rel.                     100.000000
## Eubacterium limosum et rel.                      0.000000
## Eubacterium rectale et rel.                    100.000000
## Eubacterium siraeum et rel.                      4.545455
## Eubacterium ventriosum et rel.                  93.181818
## Faecalibacterium prausnitzii et rel.           100.000000
## Fusobacteria                                     0.000000
## Gemella                                          0.000000
## Granulicatella                                   0.000000
## Haemophilus                                      0.000000
## Helicobacter                                     0.000000
## Klebisiella pneumoniae et rel.                   0.000000
## Lachnobacillus bovis et rel.                    95.454545
## Lachnospira pectinoschiza et rel.               95.454545
## Lactobacillus catenaformis et rel.               0.000000
## Lactobacillus gasseri et rel.                    0.000000
## Lactobacillus plantarum et rel.                  0.000000
## Lactobacillus salivarius et rel.                 0.000000
## Lactococcus                                     11.363636
## Leminorella                                      0.000000
## Megamonas hypermegale et rel.                    0.000000
## Megasphaera elsdenii et rel.                     0.000000
## Methylobacterium                                 0.000000
## Micrococcaceae                                   0.000000
## Mitsuokella multiacida et rel.                   2.272727
## Moraxellaceae                                    0.000000
## Novosphingobium                                  0.000000
## Oceanospirillum                                  0.000000
## Oscillospira guillermondii et rel.              81.818182
## Outgrouping clostridium cluster XIVa            88.636364
## Oxalobacter formigenes et rel.                  43.181818
## Papillibacter cinnamivorans et rel.             95.454545
## Parabacteroides distasonis et rel.              52.272727
## Peptococcus niger et rel.                        0.000000
## Peptostreptococcus anaerobius et rel.            0.000000
## Peptostreptococcus micros et rel.                0.000000
## Phascolarctobacterium faecium et rel.            0.000000
## Prevotella melaninogenica et rel.               13.636364
## Prevotella oralis et rel.                        9.090909
## Prevotella ruminicola et rel.                    2.272727
## Prevotella tannerae et rel.                     38.636364
## Propionibacterium                                0.000000
## Proteus et rel.                                  0.000000
## Pseudomonas                                      2.272727
## Roseburia intestinalis et rel.                 100.000000
## Ruminococcus bromii et rel.                     79.545455
## Ruminococcus callidus et rel.                   38.636364
## Ruminococcus gnavus et rel.                     97.727273
## Ruminococcus lactaris et rel.                  100.000000
## Ruminococcus obeum et rel.                     100.000000
## Serratia                                         4.545455
## Sporobacter termitidis et rel.                  36.363636
## Staphylococcus                                   0.000000
## Streptococcus bovis et rel.                     95.454545
## Streptococcus intermedius et rel.               25.000000
## Streptococcus mitis et rel.                     90.909091
## Subdoligranulum variable at rel.                93.181818
## Sutterella wadsworthia et rel.                   0.000000
## Tannerella et rel.                               0.000000
## Uncultured Bacteroidetes                        15.909091
## Uncultured Chroococcales                         2.272727
## Uncultured Clostridiales I                      18.181818
## Uncultured Clostridiales II                      0.000000
## Uncultured Mollicutes                            2.272727
## Uncultured Selenomonadaceae                      0.000000
## Veillonella                                      0.000000
## Weissella et rel.                                0.000000
## Vibrio                                           0.000000
## Wissella et rel.                                 0.000000
## Xanthomonadaceae                                 0.000000
## Yersinia et rel.                                 0.000000
##                                       0.00799626067236486
## Actinomycetaceae                                 0.000000
## Aerococcus                                       0.000000
## Aeromonas                                        0.000000
## Akkermansia                                     70.454545
## Alcaligenes faecalis et rel.                     0.000000
## Allistipes et rel.                              54.545455
## Anaerobiospirillum                               0.000000
## Anaerofustis                                     0.000000
## Anaerostipes caccae et rel.                    100.000000
## Anaerotruncus colihominis et rel.                0.000000
## Anaerovorax odorimutans et rel.                  2.272727
## Aneurinibacillus                                 0.000000
## Aquabacterium                                    2.272727
## Asteroleplasma et rel.                           0.000000
## Atopobium                                        0.000000
## Bacillus                                         0.000000
## Bacteroides fragilis et rel.                    38.636364
## Bacteroides intestinalis et rel.                70.454545
## Bacteroides ovatus et rel.                      15.909091
## Bacteroides plebeius et rel.                    15.909091
## Bacteroides splachnicus et rel.                  2.272727
## Bacteroides stercoris et rel.                   18.181818
## Bacteroides uniformis et rel.                   72.727273
## Bacteroides vulgatus et rel.                    88.636364
## Bifidobacterium                                 88.636364
## Bilophila et rel.                                0.000000
## Brachyspira                                      0.000000
## Bryantella formatexigens et rel.                90.909091
## Bulleidia moorei et rel.                         0.000000
## Burkholderia                                     0.000000
## Butyrivibrio crossotus et rel.                  25.000000
## Campylobacter                                    0.000000
## Catenibacterium mitsuokai et rel.                2.272727
## Clostridium (sensu stricto)                      0.000000
## Clostridium cellulosi et rel.                   34.090909
## Clostridium colinum et rel.                      2.272727
## Clostridium difficile et rel.                    0.000000
## Clostridium felsineum et rel.                    0.000000
## Clostridium leptum et rel.                      56.818182
## Clostridium nexile et rel.                     100.000000
## Clostridium orbiscindens et rel.                56.818182
## Clostridium ramosum et rel.                      0.000000
## Clostridium sphenoides et rel.                  90.909091
## Clostridium stercorarium et rel.                 0.000000
## Clostridium symbiosum et rel.                   95.454545
## Clostridium thermocellum et rel.                 0.000000
## Collinsella                                      4.545455
## Coprobacillus catenaformis et rel.               0.000000
## Coprococcus eutactus et rel.                    97.727273
## Corynebacterium                                  0.000000
## Desulfovibrio et rel.                            0.000000
## Dialister                                       13.636364
## Dorea formicigenerans et rel.                  100.000000
## Eggerthella lenta et rel.                        0.000000
## Enterobacter aerogenes et rel.                   0.000000
## Enterococcus                                     0.000000
## Escherichia coli et rel.                         2.272727
## Eubacterium biforme et rel.                     13.636364
## Eubacterium cylindroides et rel.                 2.272727
## Eubacterium hallii et rel.                     100.000000
## Eubacterium limosum et rel.                      0.000000
## Eubacterium rectale et rel.                     93.181818
## Eubacterium siraeum et rel.                      0.000000
## Eubacterium ventriosum et rel.                  88.636364
## Faecalibacterium prausnitzii et rel.           100.000000
## Fusobacteria                                     0.000000
## Gemella                                          0.000000
## Granulicatella                                   0.000000
## Haemophilus                                      0.000000
## Helicobacter                                     0.000000
## Klebisiella pneumoniae et rel.                   0.000000
## Lachnobacillus bovis et rel.                    84.090909
## Lachnospira pectinoschiza et rel.               84.090909
## Lactobacillus catenaformis et rel.               0.000000
## Lactobacillus gasseri et rel.                    0.000000
## Lactobacillus plantarum et rel.                  0.000000
## Lactobacillus salivarius et rel.                 0.000000
## Lactococcus                                      6.818182
## Leminorella                                      0.000000
## Megamonas hypermegale et rel.                    0.000000
## Megasphaera elsdenii et rel.                     0.000000
## Methylobacterium                                 0.000000
## Micrococcaceae                                   0.000000
## Mitsuokella multiacida et rel.                   2.272727
## Moraxellaceae                                    0.000000
## Novosphingobium                                  0.000000
## Oceanospirillum                                  0.000000
## Oscillospira guillermondii et rel.              56.818182
## Outgrouping clostridium cluster XIVa            65.909091
## Oxalobacter formigenes et rel.                  34.090909
## Papillibacter cinnamivorans et rel.             90.909091
## Parabacteroides distasonis et rel.              36.363636
## Peptococcus niger et rel.                        0.000000
## Peptostreptococcus anaerobius et rel.            0.000000
## Peptostreptococcus micros et rel.                0.000000
## Phascolarctobacterium faecium et rel.            0.000000
## Prevotella melaninogenica et rel.               13.636364
## Prevotella oralis et rel.                        2.272727
## Prevotella ruminicola et rel.                    0.000000
## Prevotella tannerae et rel.                     25.000000
## Propionibacterium                                0.000000
## Proteus et rel.                                  0.000000
## Pseudomonas                                      2.272727
## Roseburia intestinalis et rel.                 100.000000
## Ruminococcus bromii et rel.                     77.272727
## Ruminococcus callidus et rel.                   25.000000
## Ruminococcus gnavus et rel.                     79.545455
## Ruminococcus lactaris et rel.                   97.727273
## Ruminococcus obeum et rel.                      97.727273
## Serratia                                         2.272727
## Sporobacter termitidis et rel.                  29.545455
## Staphylococcus                                   0.000000
## Streptococcus bovis et rel.                     84.090909
## Streptococcus intermedius et rel.                6.818182
## Streptococcus mitis et rel.                     72.727273
## Subdoligranulum variable at rel.                93.181818
## Sutterella wadsworthia et rel.                   0.000000
## Tannerella et rel.                               0.000000
## Uncultured Bacteroidetes                        11.363636
## Uncultured Chroococcales                         2.272727
## Uncultured Clostridiales I                      13.636364
## Uncultured Clostridiales II                      0.000000
## Uncultured Mollicutes                            0.000000
## Uncultured Selenomonadaceae                      0.000000
## Veillonella                                      0.000000
## Weissella et rel.                                0.000000
## Vibrio                                           0.000000
## Wissella et rel.                                 0.000000
## Xanthomonadaceae                                 0.000000
## Yersinia et rel.                                 0.000000
##                                       0.0107614500053629
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                    47.727273
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                             29.545455
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                    86.363636
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 2.272727
## Aneurinibacillus                                0.000000
## Aquabacterium                                   2.272727
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                   27.272727
## Bacteroides intestinalis et rel.               52.272727
## Bacteroides ovatus et rel.                      6.818182
## Bacteroides plebeius et rel.                    6.818182
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                  13.636364
## Bacteroides uniformis et rel.                  70.454545
## Bacteroides vulgatus et rel.                   81.818182
## Bifidobacterium                                84.090909
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.               70.454545
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  6.818182
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               2.272727
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                  20.454545
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                     43.181818
## Clostridium nexile et rel.                    100.000000
## Clostridium orbiscindens et rel.               25.000000
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                 84.090909
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                  86.363636
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                   97.727273
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       4.545455
## Dorea formicigenerans et rel.                  93.181818
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     9.090909
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                    100.000000
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                    86.363636
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                 81.818182
## Faecalibacterium prausnitzii et rel.           93.181818
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                   70.454545
## Lachnospira pectinoschiza et rel.              75.000000
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     2.272727
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.             43.181818
## Outgrouping clostridium cluster XIVa           43.181818
## Oxalobacter formigenes et rel.                 22.727273
## Papillibacter cinnamivorans et rel.            84.090909
## Parabacteroides distasonis et rel.             11.363636
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.              13.636364
## Prevotella oralis et rel.                       2.272727
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                    13.636364
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 97.727273
## Ruminococcus bromii et rel.                    61.363636
## Ruminococcus callidus et rel.                  13.636364
## Ruminococcus gnavus et rel.                    63.636364
## Ruminococcus lactaris et rel.                  90.909091
## Ruminococcus obeum et rel.                     95.454545
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                 22.727273
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                    63.636364
## Streptococcus intermedius et rel.               4.545455
## Streptococcus mitis et rel.                    59.090909
## Subdoligranulum variable at rel.               84.090909
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        9.090909
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      2.272727
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0144828703018852
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                    36.363636
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                              9.090909
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                    75.000000
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 0.000000
## Aneurinibacillus                                0.000000
## Aquabacterium                                   0.000000
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                   11.363636
## Bacteroides intestinalis et rel.               34.090909
## Bacteroides ovatus et rel.                      2.272727
## Bacteroides plebeius et rel.                    6.818182
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                   6.818182
## Bacteroides uniformis et rel.                  61.363636
## Bacteroides vulgatus et rel.                   81.818182
## Bifidobacterium                                70.454545
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.               40.909091
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  0.000000
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               2.272727
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                   6.818182
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                     31.818182
## Clostridium nexile et rel.                     93.181818
## Clostridium orbiscindens et rel.                4.545455
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                 65.909091
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                  70.454545
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                   95.454545
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       0.000000
## Dorea formicigenerans et rel.                  84.090909
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     4.545455
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                    100.000000
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                    70.454545
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                 59.090909
## Faecalibacterium prausnitzii et rel.           88.636364
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                   29.545455
## Lachnospira pectinoschiza et rel.              38.636364
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     2.272727
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.             34.090909
## Outgrouping clostridium cluster XIVa           27.272727
## Oxalobacter formigenes et rel.                  9.090909
## Papillibacter cinnamivorans et rel.            68.181818
## Parabacteroides distasonis et rel.              4.545455
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.              11.363636
## Prevotella oralis et rel.                       0.000000
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                     4.545455
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 88.636364
## Ruminococcus bromii et rel.                    50.000000
## Ruminococcus callidus et rel.                   4.545455
## Ruminococcus gnavus et rel.                    27.272727
## Ruminococcus lactaris et rel.                  68.181818
## Ruminococcus obeum et rel.                     95.454545
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                  6.818182
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                    47.727273
## Streptococcus intermedius et rel.               0.000000
## Streptococcus mitis et rel.                    43.181818
## Subdoligranulum variable at rel.               81.818182
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        6.818182
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      2.272727
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0194911960820056
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                    27.272727
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                              0.000000
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                    47.727273
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 0.000000
## Aneurinibacillus                                0.000000
## Aquabacterium                                   0.000000
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                    6.818182
## Bacteroides intestinalis et rel.               29.545455
## Bacteroides ovatus et rel.                      0.000000
## Bacteroides plebeius et rel.                    6.818182
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                   2.272727
## Bacteroides uniformis et rel.                  59.090909
## Bacteroides vulgatus et rel.                   77.272727
## Bifidobacterium                                65.909091
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.                9.090909
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  0.000000
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               2.272727
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                   2.272727
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                     25.000000
## Clostridium nexile et rel.                     77.272727
## Clostridium orbiscindens et rel.                0.000000
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                 38.636364
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                  31.818182
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                   88.636364
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       0.000000
## Dorea formicigenerans et rel.                  54.545455
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     0.000000
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                    100.000000
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                    36.363636
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                 40.909091
## Faecalibacterium prausnitzii et rel.           81.818182
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                   13.636364
## Lachnospira pectinoschiza et rel.              11.363636
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     2.272727
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.             25.000000
## Outgrouping clostridium cluster XIVa           13.636364
## Oxalobacter formigenes et rel.                  0.000000
## Papillibacter cinnamivorans et rel.            43.181818
## Parabacteroides distasonis et rel.              2.272727
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.              11.363636
## Prevotella oralis et rel.                       0.000000
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                     2.272727
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 79.545455
## Ruminococcus bromii et rel.                    43.181818
## Ruminococcus callidus et rel.                   0.000000
## Ruminococcus gnavus et rel.                    13.636364
## Ruminococcus lactaris et rel.                  52.272727
## Ruminococcus obeum et rel.                     95.454545
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                  0.000000
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                    31.818182
## Streptococcus intermedius et rel.               0.000000
## Streptococcus mitis et rel.                    13.636364
## Subdoligranulum variable at rel.               70.454545
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        6.818182
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      0.000000
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0262314525220694
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                    22.727273
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                              0.000000
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                    11.363636
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 0.000000
## Aneurinibacillus                                0.000000
## Aquabacterium                                   0.000000
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                    2.272727
## Bacteroides intestinalis et rel.               25.000000
## Bacteroides ovatus et rel.                      0.000000
## Bacteroides plebeius et rel.                    4.545455
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                   0.000000
## Bacteroides uniformis et rel.                  34.090909
## Bacteroides vulgatus et rel.                   70.454545
## Bifidobacterium                                59.090909
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.                2.272727
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  0.000000
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               0.000000
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                   0.000000
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                     11.363636
## Clostridium nexile et rel.                     52.272727
## Clostridium orbiscindens et rel.                0.000000
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                  9.090909
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                   9.090909
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                   81.818182
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       0.000000
## Dorea formicigenerans et rel.                  22.727273
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     0.000000
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                     97.727273
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                     6.818182
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                 25.000000
## Faecalibacterium prausnitzii et rel.           79.545455
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                    6.818182
## Lachnospira pectinoschiza et rel.               0.000000
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     2.272727
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.              9.090909
## Outgrouping clostridium cluster XIVa            4.545455
## Oxalobacter formigenes et rel.                  0.000000
## Papillibacter cinnamivorans et rel.            18.181818
## Parabacteroides distasonis et rel.              0.000000
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.              11.363636
## Prevotella oralis et rel.                       0.000000
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                     0.000000
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 54.545455
## Ruminococcus bromii et rel.                    36.363636
## Ruminococcus callidus et rel.                   0.000000
## Ruminococcus gnavus et rel.                     4.545455
## Ruminococcus lactaris et rel.                  36.363636
## Ruminococcus obeum et rel.                     95.454545
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                  0.000000
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                     9.090909
## Streptococcus intermedius et rel.               0.000000
## Streptococcus mitis et rel.                     6.818182
## Subdoligranulum variable at rel.               50.000000
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        4.545455
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      0.000000
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0353025590898871
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                    13.636364
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                              0.000000
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                     2.272727
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 0.000000
## Aneurinibacillus                                0.000000
## Aquabacterium                                   0.000000
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                    0.000000
## Bacteroides intestinalis et rel.               13.636364
## Bacteroides ovatus et rel.                      0.000000
## Bacteroides plebeius et rel.                    2.272727
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                   0.000000
## Bacteroides uniformis et rel.                  20.454545
## Bacteroides vulgatus et rel.                   63.636364
## Bifidobacterium                                40.909091
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.                0.000000
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  0.000000
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               0.000000
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                   0.000000
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                      0.000000
## Clostridium nexile et rel.                     22.727273
## Clostridium orbiscindens et rel.                0.000000
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                  0.000000
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                   2.272727
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                   47.727273
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       0.000000
## Dorea formicigenerans et rel.                   4.545455
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     0.000000
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                     86.363636
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                     2.272727
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                 15.909091
## Faecalibacterium prausnitzii et rel.           72.727273
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                    0.000000
## Lachnospira pectinoschiza et rel.               0.000000
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     0.000000
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.              0.000000
## Outgrouping clostridium cluster XIVa            4.545455
## Oxalobacter formigenes et rel.                  0.000000
## Papillibacter cinnamivorans et rel.             4.545455
## Parabacteroides distasonis et rel.              0.000000
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.               9.090909
## Prevotella oralis et rel.                       0.000000
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                     0.000000
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 40.909091
## Ruminococcus bromii et rel.                    34.090909
## Ruminococcus callidus et rel.                   0.000000
## Ruminococcus gnavus et rel.                     0.000000
## Ruminococcus lactaris et rel.                  29.545455
## Ruminococcus obeum et rel.                     88.636364
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                  0.000000
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                     4.545455
## Streptococcus intermedius et rel.               0.000000
## Streptococcus mitis et rel.                     4.545455
## Subdoligranulum variable at rel.               29.545455
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        2.272727
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      0.000000
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0475105477764315
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                    11.363636
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                              0.000000
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                     0.000000
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 0.000000
## Aneurinibacillus                                0.000000
## Aquabacterium                                   0.000000
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                    0.000000
## Bacteroides intestinalis et rel.                4.545455
## Bacteroides ovatus et rel.                      0.000000
## Bacteroides plebeius et rel.                    0.000000
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                   0.000000
## Bacteroides uniformis et rel.                  11.363636
## Bacteroides vulgatus et rel.                   59.090909
## Bifidobacterium                                25.000000
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.                0.000000
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  0.000000
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               0.000000
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                   0.000000
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                      0.000000
## Clostridium nexile et rel.                      4.545455
## Clostridium orbiscindens et rel.                0.000000
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                  0.000000
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                   0.000000
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                   15.909091
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       0.000000
## Dorea formicigenerans et rel.                   2.272727
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     0.000000
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                     68.181818
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                     0.000000
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                  0.000000
## Faecalibacterium prausnitzii et rel.           65.909091
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                    0.000000
## Lachnospira pectinoschiza et rel.               0.000000
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     0.000000
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.              0.000000
## Outgrouping clostridium cluster XIVa            0.000000
## Oxalobacter formigenes et rel.                  0.000000
## Papillibacter cinnamivorans et rel.             0.000000
## Parabacteroides distasonis et rel.              0.000000
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.               4.545455
## Prevotella oralis et rel.                       0.000000
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                     0.000000
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 20.454545
## Ruminococcus bromii et rel.                    27.272727
## Ruminococcus callidus et rel.                   0.000000
## Ruminococcus gnavus et rel.                     0.000000
## Ruminococcus lactaris et rel.                   9.090909
## Ruminococcus obeum et rel.                     84.090909
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                  0.000000
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                     0.000000
## Streptococcus intermedius et rel.               0.000000
## Streptococcus mitis et rel.                     0.000000
## Subdoligranulum variable at rel.               11.363636
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        2.272727
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      0.000000
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0639401847404089
## Actinomycetaceae                                0.000000
## Aerococcus                                      0.000000
## Aeromonas                                       0.000000
## Akkermansia                                     6.818182
## Alcaligenes faecalis et rel.                    0.000000
## Allistipes et rel.                              0.000000
## Anaerobiospirillum                              0.000000
## Anaerofustis                                    0.000000
## Anaerostipes caccae et rel.                     0.000000
## Anaerotruncus colihominis et rel.               0.000000
## Anaerovorax odorimutans et rel.                 0.000000
## Aneurinibacillus                                0.000000
## Aquabacterium                                   0.000000
## Asteroleplasma et rel.                          0.000000
## Atopobium                                       0.000000
## Bacillus                                        0.000000
## Bacteroides fragilis et rel.                    0.000000
## Bacteroides intestinalis et rel.                2.272727
## Bacteroides ovatus et rel.                      0.000000
## Bacteroides plebeius et rel.                    0.000000
## Bacteroides splachnicus et rel.                 0.000000
## Bacteroides stercoris et rel.                   0.000000
## Bacteroides uniformis et rel.                   2.272727
## Bacteroides vulgatus et rel.                   47.727273
## Bifidobacterium                                 9.090909
## Bilophila et rel.                               0.000000
## Brachyspira                                     0.000000
## Bryantella formatexigens et rel.                0.000000
## Bulleidia moorei et rel.                        0.000000
## Burkholderia                                    0.000000
## Butyrivibrio crossotus et rel.                  0.000000
## Campylobacter                                   0.000000
## Catenibacterium mitsuokai et rel.               0.000000
## Clostridium (sensu stricto)                     0.000000
## Clostridium cellulosi et rel.                   0.000000
## Clostridium colinum et rel.                     0.000000
## Clostridium difficile et rel.                   0.000000
## Clostridium felsineum et rel.                   0.000000
## Clostridium leptum et rel.                      0.000000
## Clostridium nexile et rel.                      2.272727
## Clostridium orbiscindens et rel.                0.000000
## Clostridium ramosum et rel.                     0.000000
## Clostridium sphenoides et rel.                  0.000000
## Clostridium stercorarium et rel.                0.000000
## Clostridium symbiosum et rel.                   0.000000
## Clostridium thermocellum et rel.                0.000000
## Collinsella                                     0.000000
## Coprobacillus catenaformis et rel.              0.000000
## Coprococcus eutactus et rel.                    4.545455
## Corynebacterium                                 0.000000
## Desulfovibrio et rel.                           0.000000
## Dialister                                       0.000000
## Dorea formicigenerans et rel.                   0.000000
## Eggerthella lenta et rel.                       0.000000
## Enterobacter aerogenes et rel.                  0.000000
## Enterococcus                                    0.000000
## Escherichia coli et rel.                        0.000000
## Eubacterium biforme et rel.                     0.000000
## Eubacterium cylindroides et rel.                0.000000
## Eubacterium hallii et rel.                     38.636364
## Eubacterium limosum et rel.                     0.000000
## Eubacterium rectale et rel.                     0.000000
## Eubacterium siraeum et rel.                     0.000000
## Eubacterium ventriosum et rel.                  0.000000
## Faecalibacterium prausnitzii et rel.           54.545455
## Fusobacteria                                    0.000000
## Gemella                                         0.000000
## Granulicatella                                  0.000000
## Haemophilus                                     0.000000
## Helicobacter                                    0.000000
## Klebisiella pneumoniae et rel.                  0.000000
## Lachnobacillus bovis et rel.                    0.000000
## Lachnospira pectinoschiza et rel.               0.000000
## Lactobacillus catenaformis et rel.              0.000000
## Lactobacillus gasseri et rel.                   0.000000
## Lactobacillus plantarum et rel.                 0.000000
## Lactobacillus salivarius et rel.                0.000000
## Lactococcus                                     0.000000
## Leminorella                                     0.000000
## Megamonas hypermegale et rel.                   0.000000
## Megasphaera elsdenii et rel.                    0.000000
## Methylobacterium                                0.000000
## Micrococcaceae                                  0.000000
## Mitsuokella multiacida et rel.                  2.272727
## Moraxellaceae                                   0.000000
## Novosphingobium                                 0.000000
## Oceanospirillum                                 0.000000
## Oscillospira guillermondii et rel.              0.000000
## Outgrouping clostridium cluster XIVa            0.000000
## Oxalobacter formigenes et rel.                  0.000000
## Papillibacter cinnamivorans et rel.             0.000000
## Parabacteroides distasonis et rel.              0.000000
## Peptococcus niger et rel.                       0.000000
## Peptostreptococcus anaerobius et rel.           0.000000
## Peptostreptococcus micros et rel.               0.000000
## Phascolarctobacterium faecium et rel.           0.000000
## Prevotella melaninogenica et rel.               4.545455
## Prevotella oralis et rel.                       0.000000
## Prevotella ruminicola et rel.                   0.000000
## Prevotella tannerae et rel.                     0.000000
## Propionibacterium                               0.000000
## Proteus et rel.                                 0.000000
## Pseudomonas                                     0.000000
## Roseburia intestinalis et rel.                 11.363636
## Ruminococcus bromii et rel.                    20.454545
## Ruminococcus callidus et rel.                   0.000000
## Ruminococcus gnavus et rel.                     0.000000
## Ruminococcus lactaris et rel.                   2.272727
## Ruminococcus obeum et rel.                     52.272727
## Serratia                                        0.000000
## Sporobacter termitidis et rel.                  0.000000
## Staphylococcus                                  0.000000
## Streptococcus bovis et rel.                     0.000000
## Streptococcus intermedius et rel.               0.000000
## Streptococcus mitis et rel.                     0.000000
## Subdoligranulum variable at rel.                6.818182
## Sutterella wadsworthia et rel.                  0.000000
## Tannerella et rel.                              0.000000
## Uncultured Bacteroidetes                        2.272727
## Uncultured Chroococcales                        0.000000
## Uncultured Clostridiales I                      0.000000
## Uncultured Clostridiales II                     0.000000
## Uncultured Mollicutes                           0.000000
## Uncultured Selenomonadaceae                     0.000000
## Veillonella                                     0.000000
## Weissella et rel.                               0.000000
## Vibrio                                          0.000000
## Wissella et rel.                                0.000000
## Xanthomonadaceae                                0.000000
## Yersinia et rel.                                0.000000
##                                       0.0860513594555042 0.115808806217926
## Actinomycetaceae                                0.000000          0.000000
## Aerococcus                                      0.000000          0.000000
## Aeromonas                                       0.000000          0.000000
## Akkermansia                                     2.272727          0.000000
## Alcaligenes faecalis et rel.                    0.000000          0.000000
## Allistipes et rel.                              0.000000          0.000000
## Anaerobiospirillum                              0.000000          0.000000
## Anaerofustis                                    0.000000          0.000000
## Anaerostipes caccae et rel.                     0.000000          0.000000
## Anaerotruncus colihominis et rel.               0.000000          0.000000
## Anaerovorax odorimutans et rel.                 0.000000          0.000000
## Aneurinibacillus                                0.000000          0.000000
## Aquabacterium                                   0.000000          0.000000
## Asteroleplasma et rel.                          0.000000          0.000000
## Atopobium                                       0.000000          0.000000
## Bacillus                                        0.000000          0.000000
## Bacteroides fragilis et rel.                    0.000000          0.000000
## Bacteroides intestinalis et rel.                2.272727          0.000000
## Bacteroides ovatus et rel.                      0.000000          0.000000
## Bacteroides plebeius et rel.                    0.000000          0.000000
## Bacteroides splachnicus et rel.                 0.000000          0.000000
## Bacteroides stercoris et rel.                   0.000000          0.000000
## Bacteroides uniformis et rel.                   0.000000          0.000000
## Bacteroides vulgatus et rel.                   34.090909         31.818182
## Bifidobacterium                                 4.545455          0.000000
## Bilophila et rel.                               0.000000          0.000000
## Brachyspira                                     0.000000          0.000000
## Bryantella formatexigens et rel.                0.000000          0.000000
## Bulleidia moorei et rel.                        0.000000          0.000000
## Burkholderia                                    0.000000          0.000000
## Butyrivibrio crossotus et rel.                  0.000000          0.000000
## Campylobacter                                   0.000000          0.000000
## Catenibacterium mitsuokai et rel.               0.000000          0.000000
## Clostridium (sensu stricto)                     0.000000          0.000000
## Clostridium cellulosi et rel.                   0.000000          0.000000
## Clostridium colinum et rel.                     0.000000          0.000000
## Clostridium difficile et rel.                   0.000000          0.000000
## Clostridium felsineum et rel.                   0.000000          0.000000
## Clostridium leptum et rel.                      0.000000          0.000000
## Clostridium nexile et rel.                      0.000000          0.000000
## Clostridium orbiscindens et rel.                0.000000          0.000000
## Clostridium ramosum et rel.                     0.000000          0.000000
## Clostridium sphenoides et rel.                  0.000000          0.000000
## Clostridium stercorarium et rel.                0.000000          0.000000
## Clostridium symbiosum et rel.                   0.000000          0.000000
## Clostridium thermocellum et rel.                0.000000          0.000000
## Collinsella                                     0.000000          0.000000
## Coprobacillus catenaformis et rel.              0.000000          0.000000
## Coprococcus eutactus et rel.                    0.000000          0.000000
## Corynebacterium                                 0.000000          0.000000
## Desulfovibrio et rel.                           0.000000          0.000000
## Dialister                                       0.000000          0.000000
## Dorea formicigenerans et rel.                   0.000000          0.000000
## Eggerthella lenta et rel.                       0.000000          0.000000
## Enterobacter aerogenes et rel.                  0.000000          0.000000
## Enterococcus                                    0.000000          0.000000
## Escherichia coli et rel.                        0.000000          0.000000
## Eubacterium biforme et rel.                     0.000000          0.000000
## Eubacterium cylindroides et rel.                0.000000          0.000000
## Eubacterium hallii et rel.                     15.909091          0.000000
## Eubacterium limosum et rel.                     0.000000          0.000000
## Eubacterium rectale et rel.                     0.000000          0.000000
## Eubacterium siraeum et rel.                     0.000000          0.000000
## Eubacterium ventriosum et rel.                  0.000000          0.000000
## Faecalibacterium prausnitzii et rel.           38.636364         20.454545
## Fusobacteria                                    0.000000          0.000000
## Gemella                                         0.000000          0.000000
## Granulicatella                                  0.000000          0.000000
## Haemophilus                                     0.000000          0.000000
## Helicobacter                                    0.000000          0.000000
## Klebisiella pneumoniae et rel.                  0.000000          0.000000
## Lachnobacillus bovis et rel.                    0.000000          0.000000
## Lachnospira pectinoschiza et rel.               0.000000          0.000000
## Lactobacillus catenaformis et rel.              0.000000          0.000000
## Lactobacillus gasseri et rel.                   0.000000          0.000000
## Lactobacillus plantarum et rel.                 0.000000          0.000000
## Lactobacillus salivarius et rel.                0.000000          0.000000
## Lactococcus                                     0.000000          0.000000
## Leminorella                                     0.000000          0.000000
## Megamonas hypermegale et rel.                   0.000000          0.000000
## Megasphaera elsdenii et rel.                    0.000000          0.000000
## Methylobacterium                                0.000000          0.000000
## Micrococcaceae                                  0.000000          0.000000
## Mitsuokella multiacida et rel.                  0.000000          0.000000
## Moraxellaceae                                   0.000000          0.000000
## Novosphingobium                                 0.000000          0.000000
## Oceanospirillum                                 0.000000          0.000000
## Oscillospira guillermondii et rel.              0.000000          0.000000
## Outgrouping clostridium cluster XIVa            0.000000          0.000000
## Oxalobacter formigenes et rel.                  0.000000          0.000000
## Papillibacter cinnamivorans et rel.             0.000000          0.000000
## Parabacteroides distasonis et rel.              0.000000          0.000000
## Peptococcus niger et rel.                       0.000000          0.000000
## Peptostreptococcus anaerobius et rel.           0.000000          0.000000
## Peptostreptococcus micros et rel.               0.000000          0.000000
## Phascolarctobacterium faecium et rel.           0.000000          0.000000
## Prevotella melaninogenica et rel.               4.545455          0.000000
## Prevotella oralis et rel.                       0.000000          0.000000
## Prevotella ruminicola et rel.                   0.000000          0.000000
## Prevotella tannerae et rel.                     0.000000          0.000000
## Propionibacterium                               0.000000          0.000000
## Proteus et rel.                                 0.000000          0.000000
## Pseudomonas                                     0.000000          0.000000
## Roseburia intestinalis et rel.                  0.000000          0.000000
## Ruminococcus bromii et rel.                     6.818182          2.272727
## Ruminococcus callidus et rel.                   0.000000          0.000000
## Ruminococcus gnavus et rel.                     0.000000          0.000000
## Ruminococcus lactaris et rel.                   0.000000          0.000000
## Ruminococcus obeum et rel.                     31.818182         11.363636
## Serratia                                        0.000000          0.000000
## Sporobacter termitidis et rel.                  0.000000          0.000000
## Staphylococcus                                  0.000000          0.000000
## Streptococcus bovis et rel.                     0.000000          0.000000
## Streptococcus intermedius et rel.               0.000000          0.000000
## Streptococcus mitis et rel.                     0.000000          0.000000
## Subdoligranulum variable at rel.                0.000000          0.000000
## Sutterella wadsworthia et rel.                  0.000000          0.000000
## Tannerella et rel.                              0.000000          0.000000
## Uncultured Bacteroidetes                        2.272727          0.000000
## Uncultured Chroococcales                        0.000000          0.000000
## Uncultured Clostridiales I                      0.000000          0.000000
## Uncultured Clostridiales II                     0.000000          0.000000
## Uncultured Mollicutes                           0.000000          0.000000
## Uncultured Selenomonadaceae                     0.000000          0.000000
## Veillonella                                     0.000000          0.000000
## Weissella et rel.                               0.000000          0.000000
## Vibrio                                          0.000000          0.000000
## Wissella et rel.                                0.000000          0.000000
## Xanthomonadaceae                                0.000000          0.000000
## Yersinia et rel.                                0.000000          0.000000
##                                       0.155856684687893 0.209753532181229
## Actinomycetaceae                               0.000000          0.000000
## Aerococcus                                     0.000000          0.000000
## Aeromonas                                      0.000000          0.000000
## Akkermansia                                    0.000000          0.000000
## Alcaligenes faecalis et rel.                   0.000000          0.000000
## Allistipes et rel.                             0.000000          0.000000
## Anaerobiospirillum                             0.000000          0.000000
## Anaerofustis                                   0.000000          0.000000
## Anaerostipes caccae et rel.                    0.000000          0.000000
## Anaerotruncus colihominis et rel.              0.000000          0.000000
## Anaerovorax odorimutans et rel.                0.000000          0.000000
## Aneurinibacillus                               0.000000          0.000000
## Aquabacterium                                  0.000000          0.000000
## Asteroleplasma et rel.                         0.000000          0.000000
## Atopobium                                      0.000000          0.000000
## Bacillus                                       0.000000          0.000000
## Bacteroides fragilis et rel.                   0.000000          0.000000
## Bacteroides intestinalis et rel.               0.000000          0.000000
## Bacteroides ovatus et rel.                     0.000000          0.000000
## Bacteroides plebeius et rel.                   0.000000          0.000000
## Bacteroides splachnicus et rel.                0.000000          0.000000
## Bacteroides stercoris et rel.                  0.000000          0.000000
## Bacteroides uniformis et rel.                  0.000000          0.000000
## Bacteroides vulgatus et rel.                  18.181818          2.272727
## Bifidobacterium                                0.000000          0.000000
## Bilophila et rel.                              0.000000          0.000000
## Brachyspira                                    0.000000          0.000000
## Bryantella formatexigens et rel.               0.000000          0.000000
## Bulleidia moorei et rel.                       0.000000          0.000000
## Burkholderia                                   0.000000          0.000000
## Butyrivibrio crossotus et rel.                 0.000000          0.000000
## Campylobacter                                  0.000000          0.000000
## Catenibacterium mitsuokai et rel.              0.000000          0.000000
## Clostridium (sensu stricto)                    0.000000          0.000000
## Clostridium cellulosi et rel.                  0.000000          0.000000
## Clostridium colinum et rel.                    0.000000          0.000000
## Clostridium difficile et rel.                  0.000000          0.000000
## Clostridium felsineum et rel.                  0.000000          0.000000
## Clostridium leptum et rel.                     0.000000          0.000000
## Clostridium nexile et rel.                     0.000000          0.000000
## Clostridium orbiscindens et rel.               0.000000          0.000000
## Clostridium ramosum et rel.                    0.000000          0.000000
## Clostridium sphenoides et rel.                 0.000000          0.000000
## Clostridium stercorarium et rel.               0.000000          0.000000
## Clostridium symbiosum et rel.                  0.000000          0.000000
## Clostridium thermocellum et rel.               0.000000          0.000000
## Collinsella                                    0.000000          0.000000
## Coprobacillus catenaformis et rel.             0.000000          0.000000
## Coprococcus eutactus et rel.                   0.000000          0.000000
## Corynebacterium                                0.000000          0.000000
## Desulfovibrio et rel.                          0.000000          0.000000
## Dialister                                      0.000000          0.000000
## Dorea formicigenerans et rel.                  0.000000          0.000000
## Eggerthella lenta et rel.                      0.000000          0.000000
## Enterobacter aerogenes et rel.                 0.000000          0.000000
## Enterococcus                                   0.000000          0.000000
## Escherichia coli et rel.                       0.000000          0.000000
## Eubacterium biforme et rel.                    0.000000          0.000000
## Eubacterium cylindroides et rel.               0.000000          0.000000
## Eubacterium hallii et rel.                     0.000000          0.000000
## Eubacterium limosum et rel.                    0.000000          0.000000
## Eubacterium rectale et rel.                    0.000000          0.000000
## Eubacterium siraeum et rel.                    0.000000          0.000000
## Eubacterium ventriosum et rel.                 0.000000          0.000000
## Faecalibacterium prausnitzii et rel.           2.272727          0.000000
## Fusobacteria                                   0.000000          0.000000
## Gemella                                        0.000000          0.000000
## Granulicatella                                 0.000000          0.000000
## Haemophilus                                    0.000000          0.000000
## Helicobacter                                   0.000000          0.000000
## Klebisiella pneumoniae et rel.                 0.000000          0.000000
## Lachnobacillus bovis et rel.                   0.000000          0.000000
## Lachnospira pectinoschiza et rel.              0.000000          0.000000
## Lactobacillus catenaformis et rel.             0.000000          0.000000
## Lactobacillus gasseri et rel.                  0.000000          0.000000
## Lactobacillus plantarum et rel.                0.000000          0.000000
## Lactobacillus salivarius et rel.               0.000000          0.000000
## Lactococcus                                    0.000000          0.000000
## Leminorella                                    0.000000          0.000000
## Megamonas hypermegale et rel.                  0.000000          0.000000
## Megasphaera elsdenii et rel.                   0.000000          0.000000
## Methylobacterium                               0.000000          0.000000
## Micrococcaceae                                 0.000000          0.000000
## Mitsuokella multiacida et rel.                 0.000000          0.000000
## Moraxellaceae                                  0.000000          0.000000
## Novosphingobium                                0.000000          0.000000
## Oceanospirillum                                0.000000          0.000000
## Oscillospira guillermondii et rel.             0.000000          0.000000
## Outgrouping clostridium cluster XIVa           0.000000          0.000000
## Oxalobacter formigenes et rel.                 0.000000          0.000000
## Papillibacter cinnamivorans et rel.            0.000000          0.000000
## Parabacteroides distasonis et rel.             0.000000          0.000000
## Peptococcus niger et rel.                      0.000000          0.000000
## Peptostreptococcus anaerobius et rel.          0.000000          0.000000
## Peptostreptococcus micros et rel.              0.000000          0.000000
## Phascolarctobacterium faecium et rel.          0.000000          0.000000
## Prevotella melaninogenica et rel.              0.000000          0.000000
## Prevotella oralis et rel.                      0.000000          0.000000
## Prevotella ruminicola et rel.                  0.000000          0.000000
## Prevotella tannerae et rel.                    0.000000          0.000000
## Propionibacterium                              0.000000          0.000000
## Proteus et rel.                                0.000000          0.000000
## Pseudomonas                                    0.000000          0.000000
## Roseburia intestinalis et rel.                 0.000000          0.000000
## Ruminococcus bromii et rel.                    0.000000          0.000000
## Ruminococcus callidus et rel.                  0.000000          0.000000
## Ruminococcus gnavus et rel.                    0.000000          0.000000
## Ruminococcus lactaris et rel.                  0.000000          0.000000
## Ruminococcus obeum et rel.                     0.000000          0.000000
## Serratia                                       0.000000          0.000000
## Sporobacter termitidis et rel.                 0.000000          0.000000
## Staphylococcus                                 0.000000          0.000000
## Streptococcus bovis et rel.                    0.000000          0.000000
## Streptococcus intermedius et rel.              0.000000          0.000000
## Streptococcus mitis et rel.                    0.000000          0.000000
## Subdoligranulum variable at rel.               0.000000          0.000000
## Sutterella wadsworthia et rel.                 0.000000          0.000000
## Tannerella et rel.                             0.000000          0.000000
## Uncultured Bacteroidetes                       0.000000          0.000000
## Uncultured Chroococcales                       0.000000          0.000000
## Uncultured Clostridiales I                     0.000000          0.000000
## Uncultured Clostridiales II                    0.000000          0.000000
## Uncultured Mollicutes                          0.000000          0.000000
## Uncultured Selenomonadaceae                    0.000000          0.000000
## Veillonella                                    0.000000          0.000000
## Weissella et rel.                              0.000000          0.000000
## Vibrio                                         0.000000          0.000000
## Wissella et rel.                               0.000000          0.000000
## Xanthomonadaceae                               0.000000          0.000000
## Yersinia et rel.                               0.000000          0.000000
##                                       0.282288464884301
## Actinomycetaceae                                      0
## Aerococcus                                            0
## Aeromonas                                             0
## Akkermansia                                           0
## Alcaligenes faecalis et rel.                          0
## Allistipes et rel.                                    0
## Anaerobiospirillum                                    0
## Anaerofustis                                          0
## Anaerostipes caccae et rel.                           0
## Anaerotruncus colihominis et rel.                     0
## Anaerovorax odorimutans et rel.                       0
## Aneurinibacillus                                      0
## Aquabacterium                                         0
## Asteroleplasma et rel.                                0
## Atopobium                                             0
## Bacillus                                              0
## Bacteroides fragilis et rel.                          0
## Bacteroides intestinalis et rel.                      0
## Bacteroides ovatus et rel.                            0
## Bacteroides plebeius et rel.                          0
## Bacteroides splachnicus et rel.                       0
## Bacteroides stercoris et rel.                         0
## Bacteroides uniformis et rel.                         0
## Bacteroides vulgatus et rel.                          0
## Bifidobacterium                                       0
## Bilophila et rel.                                     0
## Brachyspira                                           0
## Bryantella formatexigens et rel.                      0
## Bulleidia moorei et rel.                              0
## Burkholderia                                          0
## Butyrivibrio crossotus et rel.                        0
## Campylobacter                                         0
## Catenibacterium mitsuokai et rel.                     0
## Clostridium (sensu stricto)                           0
## Clostridium cellulosi et rel.                         0
## Clostridium colinum et rel.                           0
## Clostridium difficile et rel.                         0
## Clostridium felsineum et rel.                         0
## Clostridium leptum et rel.                            0
## Clostridium nexile et rel.                            0
## Clostridium orbiscindens et rel.                      0
## Clostridium ramosum et rel.                           0
## Clostridium sphenoides et rel.                        0
## Clostridium stercorarium et rel.                      0
## Clostridium symbiosum et rel.                         0
## Clostridium thermocellum et rel.                      0
## Collinsella                                           0
## Coprobacillus catenaformis et rel.                    0
## Coprococcus eutactus et rel.                          0
## Corynebacterium                                       0
## Desulfovibrio et rel.                                 0
## Dialister                                             0
## Dorea formicigenerans et rel.                         0
## Eggerthella lenta et rel.                             0
## Enterobacter aerogenes et rel.                        0
## Enterococcus                                          0
## Escherichia coli et rel.                              0
## Eubacterium biforme et rel.                           0
## Eubacterium cylindroides et rel.                      0
## Eubacterium hallii et rel.                            0
## Eubacterium limosum et rel.                           0
## Eubacterium rectale et rel.                           0
## Eubacterium siraeum et rel.                           0
## Eubacterium ventriosum et rel.                        0
## Faecalibacterium prausnitzii et rel.                  0
## Fusobacteria                                          0
## Gemella                                               0
## Granulicatella                                        0
## Haemophilus                                           0
## Helicobacter                                          0
## Klebisiella pneumoniae et rel.                        0
## Lachnobacillus bovis et rel.                          0
## Lachnospira pectinoschiza et rel.                     0
## Lactobacillus catenaformis et rel.                    0
## Lactobacillus gasseri et rel.                         0
## Lactobacillus plantarum et rel.                       0
## Lactobacillus salivarius et rel.                      0
## Lactococcus                                           0
## Leminorella                                           0
## Megamonas hypermegale et rel.                         0
## Megasphaera elsdenii et rel.                          0
## Methylobacterium                                      0
## Micrococcaceae                                        0
## Mitsuokella multiacida et rel.                        0
## Moraxellaceae                                         0
## Novosphingobium                                       0
## Oceanospirillum                                       0
## Oscillospira guillermondii et rel.                    0
## Outgrouping clostridium cluster XIVa                  0
## Oxalobacter formigenes et rel.                        0
## Papillibacter cinnamivorans et rel.                   0
## Parabacteroides distasonis et rel.                    0
## Peptococcus niger et rel.                             0
## Peptostreptococcus anaerobius et rel.                 0
## Peptostreptococcus micros et rel.                     0
## Phascolarctobacterium faecium et rel.                 0
## Prevotella melaninogenica et rel.                     0
## Prevotella oralis et rel.                             0
## Prevotella ruminicola et rel.                         0
## Prevotella tannerae et rel.                           0
## Propionibacterium                                     0
## Proteus et rel.                                       0
## Pseudomonas                                           0
## Roseburia intestinalis et rel.                        0
## Ruminococcus bromii et rel.                           0
## Ruminococcus callidus et rel.                         0
## Ruminococcus gnavus et rel.                           0
## Ruminococcus lactaris et rel.                         0
## Ruminococcus obeum et rel.                            0
## Serratia                                              0
## Sporobacter termitidis et rel.                        0
## Staphylococcus                                        0
## Streptococcus bovis et rel.                           0
## Streptococcus intermedius et rel.                     0
## Streptococcus mitis et rel.                           0
## Subdoligranulum variable at rel.                      0
## Sutterella wadsworthia et rel.                        0
## Tannerella et rel.                                    0
## Uncultured Bacteroidetes                              0
## Uncultured Chroococcales                              0
## Uncultured Clostridiales I                            0
## Uncultured Clostridiales II                           0
## Uncultured Mollicutes                                 0
## Uncultured Selenomonadaceae                           0
## Veillonella                                           0
## Weissella et rel.                                     0
## Vibrio                                                0
## Wissella et rel.                                      0
## Xanthomonadaceae                                      0
## Yersinia et rel.                                      0
```
