## Filtering and pruning

The external high-quality [phyloseq package](http://joey711.github.io/phyloseq/) provides a complete set of tools for data preprocessing and filtering purposes. Download example data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling:


```r
library(microbiome)
physeq <- download_microbiome("dietswap")
```

### Sample operations


```r
# Filter samples
f1 <- filterfun_sample(topp(0.1))
wh1 <- genefilter_sample(physeq, f1, A = round(0.5 * nsamples(physeq)))

# Transform sample counts
r <- transform_sample_counts(physeq, function(x) x/sum(x))

# Sample names and variables
sample_names(physeq)
```

```
##   [1] "Sample-1"   "Sample-2"   "Sample-3"   "Sample-4"   "Sample-5"  
##   [6] "Sample-6"   "Sample-7"   "Sample-8"   "Sample-9"   "Sample-10" 
##  [11] "Sample-11"  "Sample-12"  "Sample-13"  "Sample-14"  "Sample-15" 
##  [16] "Sample-16"  "Sample-17"  "Sample-18"  "Sample-19"  "Sample-20" 
##  [21] "Sample-21"  "Sample-22"  "Sample-23"  "Sample-24"  "Sample-25" 
##  [26] "Sample-26"  "Sample-27"  "Sample-28"  "Sample-29"  "Sample-30" 
##  [31] "Sample-31"  "Sample-32"  "Sample-33"  "Sample-34"  "Sample-35" 
##  [36] "Sample-36"  "Sample-37"  "Sample-38"  "Sample-39"  "Sample-40" 
##  [41] "Sample-41"  "Sample-42"  "Sample-43"  "Sample-44"  "Sample-45" 
##  [46] "Sample-46"  "Sample-47"  "Sample-48"  "Sample-49"  "Sample-50" 
##  [51] "Sample-51"  "Sample-52"  "Sample-53"  "Sample-54"  "Sample-55" 
##  [56] "Sample-56"  "Sample-57"  "Sample-58"  "Sample-59"  "Sample-60" 
##  [61] "Sample-61"  "Sample-62"  "Sample-63"  "Sample-64"  "Sample-65" 
##  [66] "Sample-66"  "Sample-67"  "Sample-68"  "Sample-69"  "Sample-70" 
##  [71] "Sample-71"  "Sample-72"  "Sample-73"  "Sample-74"  "Sample-75" 
##  [76] "Sample-76"  "Sample-77"  "Sample-78"  "Sample-79"  "Sample-80" 
##  [81] "Sample-81"  "Sample-82"  "Sample-83"  "Sample-84"  "Sample-85" 
##  [86] "Sample-86"  "Sample-87"  "Sample-88"  "Sample-89"  "Sample-90" 
##  [91] "Sample-91"  "Sample-92"  "Sample-93"  "Sample-94"  "Sample-95" 
##  [96] "Sample-96"  "Sample-97"  "Sample-98"  "Sample-99"  "Sample-100"
## [101] "Sample-101" "Sample-102" "Sample-103" "Sample-104" "Sample-105"
## [106] "Sample-106" "Sample-107" "Sample-108" "Sample-109" "Sample-110"
## [111] "Sample-111" "Sample-112" "Sample-113" "Sample-114" "Sample-115"
## [116] "Sample-116" "Sample-117" "Sample-118" "Sample-119" "Sample-120"
## [121] "Sample-121" "Sample-122" "Sample-123" "Sample-124" "Sample-125"
## [126] "Sample-126" "Sample-127" "Sample-128" "Sample-129" "Sample-130"
## [131] "Sample-131" "Sample-132" "Sample-133" "Sample-134" "Sample-135"
## [136] "Sample-136" "Sample-137" "Sample-138" "Sample-139" "Sample-140"
## [141] "Sample-141" "Sample-142" "Sample-143" "Sample-144" "Sample-145"
## [146] "Sample-146" "Sample-147" "Sample-148" "Sample-149" "Sample-150"
## [151] "Sample-151" "Sample-152" "Sample-153" "Sample-154" "Sample-155"
## [156] "Sample-156" "Sample-157" "Sample-158" "Sample-159" "Sample-160"
## [161] "Sample-161" "Sample-162" "Sample-163" "Sample-164" "Sample-165"
## [166] "Sample-166" "Sample-167" "Sample-168" "Sample-169" "Sample-170"
## [171] "Sample-171" "Sample-172" "Sample-173" "Sample-174" "Sample-175"
## [176] "Sample-176" "Sample-177" "Sample-178" "Sample-179" "Sample-180"
## [181] "Sample-181" "Sample-182" "Sample-183" "Sample-184" "Sample-185"
## [186] "Sample-186" "Sample-187" "Sample-188" "Sample-189" "Sample-190"
## [191] "Sample-191" "Sample-192" "Sample-193" "Sample-194" "Sample-195"
## [196] "Sample-196" "Sample-197" "Sample-198" "Sample-199" "Sample-200"
## [201] "Sample-201" "Sample-202" "Sample-203" "Sample-204" "Sample-205"
## [206] "Sample-206" "Sample-207" "Sample-208" "Sample-209" "Sample-210"
## [211] "Sample-211" "Sample-212" "Sample-213" "Sample-214" "Sample-215"
## [216] "Sample-216" "Sample-217" "Sample-218" "Sample-219" "Sample-220"
## [221] "Sample-221" "Sample-222"
```

```r
# Sample sums
sample_sums(physeq)
```

```
##   Sample-1   Sample-2   Sample-3   Sample-4   Sample-5   Sample-6 
##     533779    1330516    1822706     835998    1095023    1246234 
##   Sample-7   Sample-8   Sample-9  Sample-10  Sample-11  Sample-12 
##     923662    1390057     532055    1088914    1697000    1592787 
##  Sample-13  Sample-14  Sample-15  Sample-16  Sample-17  Sample-18 
##    1054096    1072161    1176416     907865    1197154     465481 
##  Sample-19  Sample-20  Sample-21  Sample-22  Sample-23  Sample-24 
##     951907    1414691     835832     981173    1325435     901774 
##  Sample-25  Sample-26  Sample-27  Sample-28  Sample-29  Sample-30 
##    1249917     403536     987229    1127083     611728     654565 
##  Sample-31  Sample-32  Sample-33  Sample-34  Sample-35  Sample-36 
##    1032341     668019     845592     482861     781082     560486 
##  Sample-37  Sample-38  Sample-39  Sample-40  Sample-41  Sample-42 
##     894666    1323356     622431     544513     484097     316319 
##  Sample-43  Sample-44  Sample-45  Sample-46  Sample-47  Sample-48 
##     657093    1010847     858987    1294276    1017225     628157 
##  Sample-49  Sample-50  Sample-51  Sample-52  Sample-53  Sample-54 
##     908766     716655     683049     693235     530990    1485286 
##  Sample-55  Sample-56  Sample-57  Sample-58  Sample-59  Sample-60 
##    1057364     112292     770586     662477     412486     904161 
##  Sample-61  Sample-62  Sample-63  Sample-64  Sample-65  Sample-66 
##     882513     525504     740420     606865     538572     559928 
##  Sample-67  Sample-68  Sample-69  Sample-70  Sample-71  Sample-72 
##     689493     515244     943809     588612     664279     655876 
##  Sample-73  Sample-74  Sample-75  Sample-76  Sample-77  Sample-78 
##     631667     811475     368610     683339    1332703     752789 
##  Sample-79  Sample-80  Sample-81  Sample-82  Sample-83  Sample-84 
##     792505     721962    1159684     815731     673390     409563 
##  Sample-85  Sample-86  Sample-87  Sample-88  Sample-89  Sample-90 
##     694105     971267     844180     951800     894490     620728 
##  Sample-91  Sample-92  Sample-93  Sample-94  Sample-95  Sample-96 
##     950427     780016     971960     888248    1046927     893922 
##  Sample-97  Sample-98  Sample-99 Sample-100 Sample-101 Sample-102 
##     861017     869160     852016     450498    1182760     566202 
## Sample-103 Sample-104 Sample-105 Sample-106 Sample-107 Sample-108 
##     982462     839092     741029    1165799    1148335    1186664 
## Sample-109 Sample-110 Sample-111 Sample-112 Sample-113 Sample-114 
##     657081     874776    1060470    1078127     966675     989513 
## Sample-115 Sample-116 Sample-117 Sample-118 Sample-119 Sample-120 
##     974252     795508     730730    1220598     825895     673907 
## Sample-121 Sample-122 Sample-123 Sample-124 Sample-125 Sample-126 
##     668405     645616    1031083     720934    1301790     689592 
## Sample-127 Sample-128 Sample-129 Sample-130 Sample-131 Sample-132 
##     880573     750609     880159     837335     447049     914824 
## Sample-133 Sample-134 Sample-135 Sample-136 Sample-137 Sample-138 
##     938551     655318     408299     848056     546182     459487 
## Sample-139 Sample-140 Sample-141 Sample-142 Sample-143 Sample-144 
##     431419     561253     380081     873960    1168249     777908 
## Sample-145 Sample-146 Sample-147 Sample-148 Sample-149 Sample-150 
##     533171     675265     355017     780174    1003943    1055571 
## Sample-151 Sample-152 Sample-153 Sample-154 Sample-155 Sample-156 
##     791031     950794     941237     728680     943288    1320093 
## Sample-157 Sample-158 Sample-159 Sample-160 Sample-161 Sample-162 
##     705123     957801     961681    1085380    1186392     781557 
## Sample-163 Sample-164 Sample-165 Sample-166 Sample-167 Sample-168 
##     936171     349836    1321791    1081265    1118383    1004835 
## Sample-169 Sample-170 Sample-171 Sample-172 Sample-173 Sample-174 
##     964274     669725     543257     442393     853555     492795 
## Sample-175 Sample-176 Sample-177 Sample-178 Sample-179 Sample-180 
##     518357     510763     854270     614309     554082    1125410 
## Sample-181 Sample-182 Sample-183 Sample-184 Sample-185 Sample-186 
##    1122391     733984     915122     419316     477731     603341 
## Sample-187 Sample-188 Sample-189 Sample-190 Sample-191 Sample-192 
##     723635     611938     811305     756951     497837     753982 
## Sample-193 Sample-194 Sample-195 Sample-196 Sample-197 Sample-198 
##     736637     475693     282841     290892     497136     637659 
## Sample-199 Sample-200 Sample-201 Sample-202 Sample-203 Sample-204 
##     363681    1193183     641105     950352     671375     771667 
## Sample-205 Sample-206 Sample-207 Sample-208 Sample-209 Sample-210 
##     422694     966553    1196410    1459662    1134974    1126670 
## Sample-211 Sample-212 Sample-213 Sample-214 Sample-215 Sample-216 
##     942820    1249272    1365180    1201674     527332     942931 
## Sample-217 Sample-218 Sample-219 Sample-220 Sample-221 Sample-222 
##    1322774     696071    1095165    1165830     856514    1356839
```

```r
get_sample(physeq, taxa_names(physeq)[5])
```

```
##   Sample-1   Sample-2   Sample-3   Sample-4   Sample-5   Sample-6 
##         90        126        188        125        107        138 
##   Sample-7   Sample-8   Sample-9  Sample-10  Sample-11  Sample-12 
##        124        115        101        438        222        258 
##  Sample-13  Sample-14  Sample-15  Sample-16  Sample-17  Sample-18 
##        133        123         85        122        162        664 
##  Sample-19  Sample-20  Sample-21  Sample-22  Sample-23  Sample-24 
##        281        745        232        126        121        146 
##  Sample-25  Sample-26  Sample-27  Sample-28  Sample-29  Sample-30 
##        115        438        407        109         86         86 
##  Sample-31  Sample-32  Sample-33  Sample-34  Sample-35  Sample-36 
##         86        120        105         90        116         83 
##  Sample-37  Sample-38  Sample-39  Sample-40  Sample-41  Sample-42 
##        136         85         88        120        110         84 
##  Sample-43  Sample-44  Sample-45  Sample-46  Sample-47  Sample-48 
##        112        129        130        184        178        129 
##  Sample-49  Sample-50  Sample-51  Sample-52  Sample-53  Sample-54 
##        151        113        199        130         92        182 
##  Sample-55  Sample-56  Sample-57  Sample-58  Sample-59  Sample-60 
##        190        159        677        132        815        379 
##  Sample-61  Sample-62  Sample-63  Sample-64  Sample-65  Sample-66 
##        112        115        580        150        140        128 
##  Sample-67  Sample-68  Sample-69  Sample-70  Sample-71  Sample-72 
##        101        179        213        127        955       8168 
##  Sample-73  Sample-74  Sample-75  Sample-76  Sample-77  Sample-78 
##        294         86        530       1380        144         98 
##  Sample-79  Sample-80  Sample-81  Sample-82  Sample-83  Sample-84 
##         90        189        255        144        433        474 
##  Sample-85  Sample-86  Sample-87  Sample-88  Sample-89  Sample-90 
##       1269        126        139        192         89        250 
##  Sample-91  Sample-92  Sample-93  Sample-94  Sample-95  Sample-96 
##        330       1423        118         98        132        277 
##  Sample-97  Sample-98  Sample-99 Sample-100 Sample-101 Sample-102 
##        797        836        287        184        155        105 
## Sample-103 Sample-104 Sample-105 Sample-106 Sample-107 Sample-108 
##        154         92       2099        134        136        143 
## Sample-109 Sample-110 Sample-111 Sample-112 Sample-113 Sample-114 
##        104        114        112        120        101        197 
## Sample-115 Sample-116 Sample-117 Sample-118 Sample-119 Sample-120 
##        103         97        101        157        138        126 
## Sample-121 Sample-122 Sample-123 Sample-124 Sample-125 Sample-126 
##        125        100        149        130        177        100 
## Sample-127 Sample-128 Sample-129 Sample-130 Sample-131 Sample-132 
##        147        108        195        132         82        126 
## Sample-133 Sample-134 Sample-135 Sample-136 Sample-137 Sample-138 
##        101         90        101        201        120         95 
## Sample-139 Sample-140 Sample-141 Sample-142 Sample-143 Sample-144 
##         94        112         92        143        138         91 
## Sample-145 Sample-146 Sample-147 Sample-148 Sample-149 Sample-150 
##        172        125        146         97        373       1039 
## Sample-151 Sample-152 Sample-153 Sample-154 Sample-155 Sample-156 
##        229       1789      11796         94        156        313 
## Sample-157 Sample-158 Sample-159 Sample-160 Sample-161 Sample-162 
##       4051       1444        131        153        233        151 
## Sample-163 Sample-164 Sample-165 Sample-166 Sample-167 Sample-168 
##        140       2008       2265        286        259        122 
## Sample-169 Sample-170 Sample-171 Sample-172 Sample-173 Sample-174 
##        159        109       2070        297        161        144 
## Sample-175 Sample-176 Sample-177 Sample-178 Sample-179 Sample-180 
##       1188        251        142        145       6333        324 
## Sample-181 Sample-182 Sample-183 Sample-184 Sample-185 Sample-186 
##        187        104         94        177       2236        277 
## Sample-187 Sample-188 Sample-189 Sample-190 Sample-191 Sample-192 
##         87        109         99         92         91        124 
## Sample-193 Sample-194 Sample-195 Sample-196 Sample-197 Sample-198 
##        126         86        123        101        127        110 
## Sample-199 Sample-200 Sample-201 Sample-202 Sample-203 Sample-204 
##        101         84        110         95        146        100 
## Sample-205 Sample-206 Sample-207 Sample-208 Sample-209 Sample-210 
##        110         93       2042        210        479       3668 
## Sample-211 Sample-212 Sample-213 Sample-214 Sample-215 Sample-216 
##       3902        161        142        125        863        552 
## Sample-217 Sample-218 Sample-219 Sample-220 Sample-221 Sample-222 
##        703        354        340        108        136        199
```


### Variable operations


```r
sample_variables(physeq)
```

```
## [1] "subject"                "gender"                
## [3] "nationality"            "group"                 
## [5] "sample"                 "timepoint"             
## [7] "timepoint.within.group" "bmi_group"
```

```r
get_variable(physeq, sample_variables(physeq)[5])[1:10]
```

```
##  [1] "Sample-1"  "Sample-2"  "Sample-3"  "Sample-4"  "Sample-5" 
##  [6] "Sample-6"  "Sample-7"  "Sample-8"  "Sample-9"  "Sample-10"
```

```r
# Assign fields to sample metadata
# sample_data(GP)$human <- ..
```

### Taxa operations


```r
# Prune taxa
ex2 <- prune_taxa(wh1, physeq)
ex3 <- prune_taxa(taxa_sums(physeq) > 0, physeq)

# Subset taxa
pseq <- subset_taxa(pseq, Phylum == "Bacteroidetes")

# Filter taxa
f <- filter_taxa(r, function(x) var(x) > 1e-05, TRUE)

# Number of taxa
ntaxa(physeq)
```

```
## [1] 130
```

```r
# Rank names
rank_names(physeq)
```

```
## [1] "Phylum" "Genus"
```

```r
# Taxa names
taxa_names(physeq)
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

```r
get_taxa_unique(physeq, "Phylum")
```

```
##  [1] "Actinobacteria"            "Bacilli"                  
##  [3] "Proteobacteria"            "Verrucomicrobia"          
##  [5] "Bacteroidetes"             "Clostridium cluster XV"   
##  [7] "Clostridium cluster XIVa"  "Clostridium cluster IV"   
##  [9] "Clostridium cluster XI"    "Asteroleplasma"           
## [11] "Spirochaetes"              "Clostridium cluster XVI"  
## [13] "Clostridium cluster XVII"  "Clostridium cluster I"    
## [15] "Clostridium cluster XVIII" "Clostridium cluster III"  
## [17] "Clostridium cluster IX"    "Fusobacteria"             
## [19] "Clostridium cluster XIII"  "Cyanobacteria"            
## [21] "Uncultured Clostridiales"  "Uncultured Mollicutes"
```

```r
get_taxa(physeq, sample_names(physeq)[5])
```

```
##                      Actinomycetaceae 
##                                    16 
##                            Aerococcus 
##                                     1 
##                             Aeromonas 
##                                     1 
##                           Akkermansia 
##                                  1331 
##          Alcaligenes faecalis et rel. 
##                                   107 
##                    Allistipes et rel. 
##                                  3087 
##                    Anaerobiospirillum 
##                                     1 
##                          Anaerofustis 
##                                     1 
##           Anaerostipes caccae et rel. 
##                                  1467 
##     Anaerotruncus colihominis et rel. 
##                                   938 
##       Anaerovorax odorimutans et rel. 
##                                   735 
##                      Aneurinibacillus 
##                                     1 
##                         Aquabacterium 
##                                     1 
##                Asteroleplasma et rel. 
##                                     1 
##                             Atopobium 
##                                    22 
##                              Bacillus 
##                                    65 
##          Bacteroides fragilis et rel. 
##                                  2088 
##      Bacteroides intestinalis et rel. 
##                                   176 
##            Bacteroides ovatus et rel. 
##                                  1800 
##          Bacteroides plebeius et rel. 
##                                  1302 
##       Bacteroides splachnicus et rel. 
##                                  1907 
##         Bacteroides stercoris et rel. 
##                                   612 
##         Bacteroides uniformis et rel. 
##                                   320 
##          Bacteroides vulgatus et rel. 
##                                 10517 
##                       Bifidobacterium 
##                                  1390 
##                     Bilophila et rel. 
##                                    70 
##                           Brachyspira 
##                                    15 
##      Bryantella formatexigens et rel. 
##                                  3397 
##              Bulleidia moorei et rel. 
##                                   126 
##                          Burkholderia 
##                                   110 
##        Butyrivibrio crossotus et rel. 
##                                 11169 
##                         Campylobacter 
##                                   253 
##     Catenibacterium mitsuokai et rel. 
##                                    32 
##           Clostridium (sensu stricto) 
##                                  1308 
##         Clostridium cellulosi et rel. 
##                                 45935 
##           Clostridium colinum et rel. 
##                                  5441 
##         Clostridium difficile et rel. 
##                                  1604 
##         Clostridium felsineum et rel. 
##                                     1 
##            Clostridium leptum et rel. 
##                                  1834 
##            Clostridium nexile et rel. 
##                                  4252 
##      Clostridium orbiscindens et rel. 
##                                 16735 
##           Clostridium ramosum et rel. 
##                                   147 
##        Clostridium sphenoides et rel. 
##                                  5923 
##      Clostridium stercorarium et rel. 
##                                   193 
##         Clostridium symbiosum et rel. 
##                                  9799 
##      Clostridium thermocellum et rel. 
##                                     1 
##                           Collinsella 
##                                  1688 
##    Coprobacillus catenaformis et rel. 
##                                   247 
##          Coprococcus eutactus et rel. 
##                                 13337 
##                       Corynebacterium 
##                                    58 
##                 Desulfovibrio et rel. 
##                                   197 
##                             Dialister 
##                                   361 
##         Dorea formicigenerans et rel. 
##                                 11937 
##             Eggerthella lenta et rel. 
##                                   254 
##        Enterobacter aerogenes et rel. 
##                                   927 
##                          Enterococcus 
##                                   247 
##              Escherichia coli et rel. 
##                                   669 
##           Eubacterium biforme et rel. 
##                                  1091 
##      Eubacterium cylindroides et rel. 
##                                   149 
##            Eubacterium hallii et rel. 
##                                  1746 
##           Eubacterium limosum et rel. 
##                                    65 
##           Eubacterium rectale et rel. 
##                                  3946 
##           Eubacterium siraeum et rel. 
##                                   171 
##        Eubacterium ventriosum et rel. 
##                                   891 
##  Faecalibacterium prausnitzii et rel. 
##                                 87375 
##                          Fusobacteria 
##                                   299 
##                               Gemella 
##                                     1 
##                        Granulicatella 
##                                     1 
##                           Haemophilus 
##                                    46 
##                          Helicobacter 
##                                   138 
##        Klebisiella pneumoniae et rel. 
##                                   397 
##          Lachnobacillus bovis et rel. 
##                                  2066 
##     Lachnospira pectinoschiza et rel. 
##                                  4238 
##    Lactobacillus catenaformis et rel. 
##                                    21 
##         Lactobacillus gasseri et rel. 
##                                   443 
##       Lactobacillus plantarum et rel. 
##                                   723 
##      Lactobacillus salivarius et rel. 
##                                    62 
##                           Lactococcus 
##                                    75 
##                           Leminorella 
##                                     1 
##         Megamonas hypermegale et rel. 
##                                    58 
##          Megasphaera elsdenii et rel. 
##                                   450 
##                      Methylobacterium 
##                                     1 
##                        Micrococcaceae 
##                                     1 
##        Mitsuokella multiacida et rel. 
##                                    82 
##                         Moraxellaceae 
##                                    47 
##                       Novosphingobium 
##                                     1 
##                       Oceanospirillum 
##                                    96 
##    Oscillospira guillermondii et rel. 
##                                 24765 
##  Outgrouping clostridium cluster XIVa 
##                                  8720 
##        Oxalobacter formigenes et rel. 
##                                   326 
##   Papillibacter cinnamivorans et rel. 
##                                  6824 
##    Parabacteroides distasonis et rel. 
##                                  2414 
##             Peptococcus niger et rel. 
##                                   126 
## Peptostreptococcus anaerobius et rel. 
##                                     1 
##     Peptostreptococcus micros et rel. 
##                                   135 
## Phascolarctobacterium faecium et rel. 
##                                  2110 
##     Prevotella melaninogenica et rel. 
##                                596658 
##             Prevotella oralis et rel. 
##                                 94935 
##         Prevotella ruminicola et rel. 
##                                   557 
##           Prevotella tannerae et rel. 
##                                   977 
##                     Propionibacterium 
##                                    63 
##                       Proteus et rel. 
##                                   255 
##                           Pseudomonas 
##                                    19 
##        Roseburia intestinalis et rel. 
##                                  1486 
##           Ruminococcus bromii et rel. 
##                                  1863 
##         Ruminococcus callidus et rel. 
##                                  2509 
##           Ruminococcus gnavus et rel. 
##                                  2957 
##         Ruminococcus lactaris et rel. 
##                                   473 
##            Ruminococcus obeum et rel. 
##                                 28848 
##                              Serratia 
##                                    25 
##        Sporobacter termitidis et rel. 
##                                  8694 
##                        Staphylococcus 
##                                    21 
##           Streptococcus bovis et rel. 
##                                  6914 
##     Streptococcus intermedius et rel. 
##                                   338 
##           Streptococcus mitis et rel. 
##                                  2720 
##      Subdoligranulum variable at rel. 
##                                 18940 
##        Sutterella wadsworthia et rel. 
##                                  3497 
##                    Tannerella et rel. 
##                                  1268 
##              Uncultured Bacteroidetes 
##                                    71 
##              Uncultured Chroococcales 
##                                     1 
##            Uncultured Clostridiales I 
##                                  1203 
##           Uncultured Clostridiales II 
##                                  1895 
##                 Uncultured Mollicutes 
##                                   541 
##           Uncultured Selenomonadaceae 
##                                     1 
##                           Veillonella 
##                                   458 
##                                Vibrio 
##                                   143 
##                     Weissella et rel. 
##                                   125 
##                      Wissella et rel. 
##                                    16 
##                      Xanthomonadaceae 
##                                   185 
##                      Yersinia et rel. 
##                                   111
```

```r
# Taxa sums
taxa_sums(physeq)
```

```
##                      Actinomycetaceae 
##                                  4289 
##                            Aerococcus 
##                                   241 
##                             Aeromonas 
##                                   330 
##                           Akkermansia 
##                               1613341 
##          Alcaligenes faecalis et rel. 
##                                 98654 
##                    Allistipes et rel. 
##                               3513027 
##                    Anaerobiospirillum 
##                                   362 
##                          Anaerofustis 
##                                   855 
##           Anaerostipes caccae et rel. 
##                                889316 
##     Anaerotruncus colihominis et rel. 
##                               1011044 
##       Anaerovorax odorimutans et rel. 
##                                196839 
##                      Aneurinibacillus 
##                                   452 
##                         Aquabacterium 
##                                  4100 
##                Asteroleplasma et rel. 
##                                   222 
##                             Atopobium 
##                                  5171 
##                              Bacillus 
##                                 14538 
##          Bacteroides fragilis et rel. 
##                               2539567 
##      Bacteroides intestinalis et rel. 
##                                199684 
##            Bacteroides ovatus et rel. 
##                               1516522 
##          Bacteroides plebeius et rel. 
##                                596972 
##       Bacteroides splachnicus et rel. 
##                                833871 
##         Bacteroides stercoris et rel. 
##                                367524 
##         Bacteroides uniformis et rel. 
##                               1302113 
##          Bacteroides vulgatus et rel. 
##                              17613054 
##                       Bifidobacterium 
##                               1019045 
##                     Bilophila et rel. 
##                                 21142 
##                           Brachyspira 
##                                  4846 
##      Bryantella formatexigens et rel. 
##                               1238787 
##              Bulleidia moorei et rel. 
##                                 60824 
##                          Burkholderia 
##                                 13795 
##        Butyrivibrio crossotus et rel. 
##                               2556998 
##                         Campylobacter 
##                                 55937 
##     Catenibacterium mitsuokai et rel. 
##                                 60582 
##           Clostridium (sensu stricto) 
##                                352244 
##         Clostridium cellulosi et rel. 
##                               9668002 
##           Clostridium colinum et rel. 
##                                405221 
##         Clostridium difficile et rel. 
##                                400561 
##         Clostridium felsineum et rel. 
##                                   227 
##            Clostridium leptum et rel. 
##                               1162473 
##            Clostridium nexile et rel. 
##                                695112 
##      Clostridium orbiscindens et rel. 
##                               3410939 
##           Clostridium ramosum et rel. 
##                                 33691 
##        Clostridium sphenoides et rel. 
##                               1130553 
##      Clostridium stercorarium et rel. 
##                                126228 
##         Clostridium symbiosum et rel. 
##                               3614491 
##      Clostridium thermocellum et rel. 
##                                   222 
##                           Collinsella 
##                                173988 
##    Coprobacillus catenaformis et rel. 
##                                 44852 
##          Coprococcus eutactus et rel. 
##                               1218564 
##                       Corynebacterium 
##                                 13410 
##                 Desulfovibrio et rel. 
##                                 63536 
##                             Dialister 
##                               1248524 
##         Dorea formicigenerans et rel. 
##                               1445868 
##             Eggerthella lenta et rel. 
##                                 63978 
##        Enterobacter aerogenes et rel. 
##                                281964 
##                          Enterococcus 
##                                139975 
##              Escherichia coli et rel. 
##                               1122252 
##           Eubacterium biforme et rel. 
##                                380380 
##      Eubacterium cylindroides et rel. 
##                                 35333 
##            Eubacterium hallii et rel. 
##                                286840 
##           Eubacterium limosum et rel. 
##                                 16368 
##           Eubacterium rectale et rel. 
##                                726610 
##           Eubacterium siraeum et rel. 
##                                 51100 
##        Eubacterium ventriosum et rel. 
##                                278287 
##  Faecalibacterium prausnitzii et rel. 
##                               6845609 
##                          Fusobacteria 
##                                 87933 
##                               Gemella 
##                                   840 
##                        Granulicatella 
##                                  2803 
##                           Haemophilus 
##                                 33623 
##                          Helicobacter 
##                                 30274 
##        Klebisiella pneumoniae et rel. 
##                                186187 
##          Lachnobacillus bovis et rel. 
##                                728925 
##     Lachnospira pectinoschiza et rel. 
##                                897134 
##    Lactobacillus catenaformis et rel. 
##                                 11698 
##         Lactobacillus gasseri et rel. 
##                                111923 
##       Lactobacillus plantarum et rel. 
##                                260289 
##      Lactobacillus salivarius et rel. 
##                                 33920 
##                           Lactococcus 
##                                 22555 
##                           Leminorella 
##                                  6884 
##         Megamonas hypermegale et rel. 
##                                 13236 
##          Megasphaera elsdenii et rel. 
##                                599451 
##                      Methylobacterium 
##                                   222 
##                        Micrococcaceae 
##                                   222 
##        Mitsuokella multiacida et rel. 
##                                709390 
##                         Moraxellaceae 
##                                  7624 
##                       Novosphingobium 
##                                   402 
##                       Oceanospirillum 
##                                 41540 
##    Oscillospira guillermondii et rel. 
##                              22208228 
##  Outgrouping clostridium cluster XIVa 
##                                709171 
##        Oxalobacter formigenes et rel. 
##                                484519 
##   Papillibacter cinnamivorans et rel. 
##                                431703 
##    Parabacteroides distasonis et rel. 
##                               1509958 
##             Peptococcus niger et rel. 
##                                 54555 
## Peptostreptococcus anaerobius et rel. 
##                                   330 
##     Peptostreptococcus micros et rel. 
##                                 32783 
## Phascolarctobacterium faecium et rel. 
##                                460846 
##     Prevotella melaninogenica et rel. 
##                              55102081 
##             Prevotella oralis et rel. 
##                               8405109 
##         Prevotella ruminicola et rel. 
##                                 46918 
##           Prevotella tannerae et rel. 
##                                744898 
##                     Propionibacterium 
##                                 16504 
##                       Proteus et rel. 
##                                 65761 
##                           Pseudomonas 
##                                  3995 
##        Roseburia intestinalis et rel. 
##                                293290 
##           Ruminococcus bromii et rel. 
##                                750430 
##         Ruminococcus callidus et rel. 
##                                857217 
##           Ruminococcus gnavus et rel. 
##                                406358 
##         Ruminococcus lactaris et rel. 
##                                137526 
##            Ruminococcus obeum et rel. 
##                               3296629 
##                              Serratia 
##                                 62967 
##        Sporobacter termitidis et rel. 
##                               4107609 
##                        Staphylococcus 
##                                  4820 
##           Streptococcus bovis et rel. 
##                               1069640 
##     Streptococcus intermedius et rel. 
##                                 90656 
##           Streptococcus mitis et rel. 
##                                602624 
##      Subdoligranulum variable at rel. 
##                               3303469 
##        Sutterella wadsworthia et rel. 
##                                576410 
##                    Tannerella et rel. 
##                                465584 
##              Uncultured Bacteroidetes 
##                                 44595 
##              Uncultured Chroococcales 
##                                   905 
##            Uncultured Clostridiales I 
##                               1227178 
##           Uncultured Clostridiales II 
##                               1066574 
##                 Uncultured Mollicutes 
##                                628130 
##           Uncultured Selenomonadaceae 
##                                 23037 
##                           Veillonella 
##                                113879 
##                                Vibrio 
##                                 37200 
##                     Weissella et rel. 
##                                 49089 
##                      Wissella et rel. 
##                                  5052 
##                      Xanthomonadaceae 
##                                 18677 
##                      Yersinia et rel. 
##                                 33733
```


### Merging phyloseq operations


```r
# merge_phyloseq()
# merge_taxa()
# tax_glom()
```


### Rarification


```r
physeq.rarified <- rarefy_even_depth(physeq)
```

