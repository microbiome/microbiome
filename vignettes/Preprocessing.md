## Filtering and pruning

The external high-quality [phyloseq package](http://joey711.github.io/phyloseq/) provides a complete set of tools for data preprocessing and filtering purposes. Download example data from [O'Keefe et al. Nat. Comm. 6:6342, 2015](http://dx.doi.org/10.1038/ncomms7342) from [Data Dryad](http://dx.doi.org/10.5061/dryad.1mn1n). This is a two-week diet swap study between western (USA) and traditional (rural Africa) diets including microbiota profiling:


```r
library(microbiome)
pseq <- download_microbiome("dietswap")
```

### Sample operations


Transform sample counts


```r
r <- transform_sample_counts(pseq, function(x) x/sum(x))
```

Sample names and variables


```r
sample_names(pseq)
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

Sample sums


```r
sample_sums(pseq)
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

Abundance for species ‘i’ in each sample


```r
get_sample(pseq, taxa_names(pseq)[5])
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
head(get_variable(pseq, sample_variables(pseq)[1]))
```

```
## [1] byn nms olt pku qjy riv
## 38 Levels: azh azl byn byu cxj dwc dwk eve fua fud gtd gty hsf irh ... zaq
```

```r
# Assign fields to sample metadata
# sample_data(GP)$human <- ..
```

### Taxa operations


Filter samples


```r
f1 <- filterfun_sample(topp(0.1))
taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))
taxa
```

```
##                      Actinomycetaceae 
##                                 FALSE 
##                            Aerococcus 
##                                 FALSE 
##                             Aeromonas 
##                                 FALSE 
##                           Akkermansia 
##                                 FALSE 
##          Alcaligenes faecalis et rel. 
##                                 FALSE 
##                    Allistipes et rel. 
##                                 FALSE 
##                    Anaerobiospirillum 
##                                 FALSE 
##                          Anaerofustis 
##                                 FALSE 
##           Anaerostipes caccae et rel. 
##                                 FALSE 
##     Anaerotruncus colihominis et rel. 
##                                 FALSE 
##       Anaerovorax odorimutans et rel. 
##                                 FALSE 
##                      Aneurinibacillus 
##                                 FALSE 
##                         Aquabacterium 
##                                 FALSE 
##                Asteroleplasma et rel. 
##                                 FALSE 
##                             Atopobium 
##                                 FALSE 
##                              Bacillus 
##                                 FALSE 
##          Bacteroides fragilis et rel. 
##                                 FALSE 
##      Bacteroides intestinalis et rel. 
##                                 FALSE 
##            Bacteroides ovatus et rel. 
##                                 FALSE 
##          Bacteroides plebeius et rel. 
##                                 FALSE 
##       Bacteroides splachnicus et rel. 
##                                 FALSE 
##         Bacteroides stercoris et rel. 
##                                 FALSE 
##         Bacteroides uniformis et rel. 
##                                 FALSE 
##          Bacteroides vulgatus et rel. 
##                                  TRUE 
##                       Bifidobacterium 
##                                 FALSE 
##                     Bilophila et rel. 
##                                 FALSE 
##                           Brachyspira 
##                                 FALSE 
##      Bryantella formatexigens et rel. 
##                                 FALSE 
##              Bulleidia moorei et rel. 
##                                 FALSE 
##                          Burkholderia 
##                                 FALSE 
##        Butyrivibrio crossotus et rel. 
##                                 FALSE 
##                         Campylobacter 
##                                 FALSE 
##     Catenibacterium mitsuokai et rel. 
##                                 FALSE 
##           Clostridium (sensu stricto) 
##                                 FALSE 
##         Clostridium cellulosi et rel. 
##                                  TRUE 
##           Clostridium colinum et rel. 
##                                 FALSE 
##         Clostridium difficile et rel. 
##                                 FALSE 
##         Clostridium felsineum et rel. 
##                                 FALSE 
##            Clostridium leptum et rel. 
##                                 FALSE 
##            Clostridium nexile et rel. 
##                                 FALSE 
##      Clostridium orbiscindens et rel. 
##                                  TRUE 
##           Clostridium ramosum et rel. 
##                                 FALSE 
##        Clostridium sphenoides et rel. 
##                                 FALSE 
##      Clostridium stercorarium et rel. 
##                                 FALSE 
##         Clostridium symbiosum et rel. 
##                                  TRUE 
##      Clostridium thermocellum et rel. 
##                                 FALSE 
##                           Collinsella 
##                                 FALSE 
##    Coprobacillus catenaformis et rel. 
##                                 FALSE 
##          Coprococcus eutactus et rel. 
##                                 FALSE 
##                       Corynebacterium 
##                                 FALSE 
##                 Desulfovibrio et rel. 
##                                 FALSE 
##                             Dialister 
##                                 FALSE 
##         Dorea formicigenerans et rel. 
##                                 FALSE 
##             Eggerthella lenta et rel. 
##                                 FALSE 
##        Enterobacter aerogenes et rel. 
##                                 FALSE 
##                          Enterococcus 
##                                 FALSE 
##              Escherichia coli et rel. 
##                                 FALSE 
##           Eubacterium biforme et rel. 
##                                 FALSE 
##      Eubacterium cylindroides et rel. 
##                                 FALSE 
##            Eubacterium hallii et rel. 
##                                 FALSE 
##           Eubacterium limosum et rel. 
##                                 FALSE 
##           Eubacterium rectale et rel. 
##                                 FALSE 
##           Eubacterium siraeum et rel. 
##                                 FALSE 
##        Eubacterium ventriosum et rel. 
##                                 FALSE 
##  Faecalibacterium prausnitzii et rel. 
##                                  TRUE 
##                          Fusobacteria 
##                                 FALSE 
##                               Gemella 
##                                 FALSE 
##                        Granulicatella 
##                                 FALSE 
##                           Haemophilus 
##                                 FALSE 
##                          Helicobacter 
##                                 FALSE 
##        Klebisiella pneumoniae et rel. 
##                                 FALSE 
##          Lachnobacillus bovis et rel. 
##                                 FALSE 
##     Lachnospira pectinoschiza et rel. 
##                                 FALSE 
##    Lactobacillus catenaformis et rel. 
##                                 FALSE 
##         Lactobacillus gasseri et rel. 
##                                 FALSE 
##       Lactobacillus plantarum et rel. 
##                                 FALSE 
##      Lactobacillus salivarius et rel. 
##                                 FALSE 
##                           Lactococcus 
##                                 FALSE 
##                           Leminorella 
##                                 FALSE 
##         Megamonas hypermegale et rel. 
##                                 FALSE 
##          Megasphaera elsdenii et rel. 
##                                 FALSE 
##                      Methylobacterium 
##                                 FALSE 
##                        Micrococcaceae 
##                                 FALSE 
##        Mitsuokella multiacida et rel. 
##                                 FALSE 
##                         Moraxellaceae 
##                                 FALSE 
##                       Novosphingobium 
##                                 FALSE 
##                       Oceanospirillum 
##                                 FALSE 
##    Oscillospira guillermondii et rel. 
##                                  TRUE 
##  Outgrouping clostridium cluster XIVa 
##                                 FALSE 
##        Oxalobacter formigenes et rel. 
##                                 FALSE 
##   Papillibacter cinnamivorans et rel. 
##                                 FALSE 
##    Parabacteroides distasonis et rel. 
##                                 FALSE 
##             Peptococcus niger et rel. 
##                                 FALSE 
## Peptostreptococcus anaerobius et rel. 
##                                 FALSE 
##     Peptostreptococcus micros et rel. 
##                                 FALSE 
## Phascolarctobacterium faecium et rel. 
##                                 FALSE 
##     Prevotella melaninogenica et rel. 
##                                  TRUE 
##             Prevotella oralis et rel. 
##                                  TRUE 
##         Prevotella ruminicola et rel. 
##                                 FALSE 
##           Prevotella tannerae et rel. 
##                                 FALSE 
##                     Propionibacterium 
##                                 FALSE 
##                       Proteus et rel. 
##                                 FALSE 
##                           Pseudomonas 
##                                 FALSE 
##        Roseburia intestinalis et rel. 
##                                 FALSE 
##           Ruminococcus bromii et rel. 
##                                 FALSE 
##         Ruminococcus callidus et rel. 
##                                 FALSE 
##           Ruminococcus gnavus et rel. 
##                                 FALSE 
##         Ruminococcus lactaris et rel. 
##                                 FALSE 
##            Ruminococcus obeum et rel. 
##                                  TRUE 
##                              Serratia 
##                                 FALSE 
##        Sporobacter termitidis et rel. 
##                                  TRUE 
##                        Staphylococcus 
##                                 FALSE 
##           Streptococcus bovis et rel. 
##                                 FALSE 
##     Streptococcus intermedius et rel. 
##                                 FALSE 
##           Streptococcus mitis et rel. 
##                                 FALSE 
##      Subdoligranulum variable at rel. 
##                                  TRUE 
##        Sutterella wadsworthia et rel. 
##                                 FALSE 
##                    Tannerella et rel. 
##                                 FALSE 
##              Uncultured Bacteroidetes 
##                                 FALSE 
##              Uncultured Chroococcales 
##                                 FALSE 
##            Uncultured Clostridiales I 
##                                 FALSE 
##           Uncultured Clostridiales II 
##                                 FALSE 
##                 Uncultured Mollicutes 
##                                 FALSE 
##           Uncultured Selenomonadaceae 
##                                 FALSE 
##                           Veillonella 
##                                 FALSE 
##                                Vibrio 
##                                 FALSE 
##                     Weissella et rel. 
##                                 FALSE 
##                      Wissella et rel. 
##                                 FALSE 
##                      Xanthomonadaceae 
##                                 FALSE 
##                      Yersinia et rel. 
##                                 FALSE
```

Prune taxa


```r
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

Number of taxa


```r
ntaxa(pseq)
```

```
## [1] 16
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
##  [1] "Allistipes et rel."                
##  [2] "Bacteroides fragilis et rel."      
##  [3] "Bacteroides intestinalis et rel."  
##  [4] "Bacteroides ovatus et rel."        
##  [5] "Bacteroides plebeius et rel."      
##  [6] "Bacteroides splachnicus et rel."   
##  [7] "Bacteroides stercoris et rel."     
##  [8] "Bacteroides uniformis et rel."     
##  [9] "Bacteroides vulgatus et rel."      
## [10] "Parabacteroides distasonis et rel."
## [11] "Prevotella melaninogenica et rel." 
## [12] "Prevotella oralis et rel."         
## [13] "Prevotella ruminicola et rel."     
## [14] "Prevotella tannerae et rel."       
## [15] "Tannerella et rel."                
## [16] "Uncultured Bacteroidetes"
```


Pick taxa


```r
# Unique phyla
head(get_taxa_unique(pseq, "Phylum"))
```

```
## [1] "Bacteroidetes"
```

```r
# Taxa by sample 
get_taxa(pseq, sample_names(pseq)[1])
```

```
##                 Allistipes et rel.       Bacteroides fragilis et rel. 
##                              21222                              27925 
##   Bacteroides intestinalis et rel.         Bacteroides ovatus et rel. 
##                                758                              26913 
##       Bacteroides plebeius et rel.    Bacteroides splachnicus et rel. 
##                               4073                               2835 
##      Bacteroides stercoris et rel.      Bacteroides uniformis et rel. 
##                               1347                              19582 
##       Bacteroides vulgatus et rel. Parabacteroides distasonis et rel. 
##                             175035                              11772 
##  Prevotella melaninogenica et rel.          Prevotella oralis et rel. 
##                               1563                               1386 
##      Prevotella ruminicola et rel.        Prevotella tannerae et rel. 
##                                 70                               3964 
##                 Tannerella et rel.           Uncultured Bacteroidetes 
##                               3448                                113
```


Taxa sums


```r
taxa_sums(pseq)
```

```
##                 Allistipes et rel.       Bacteroides fragilis et rel. 
##                            3513027                            2539567 
##   Bacteroides intestinalis et rel.         Bacteroides ovatus et rel. 
##                             199684                            1516522 
##       Bacteroides plebeius et rel.    Bacteroides splachnicus et rel. 
##                             596972                             833871 
##      Bacteroides stercoris et rel.      Bacteroides uniformis et rel. 
##                             367524                            1302113 
##       Bacteroides vulgatus et rel. Parabacteroides distasonis et rel. 
##                           17613054                            1509958 
##  Prevotella melaninogenica et rel.          Prevotella oralis et rel. 
##                           55102081                            8405109 
##      Prevotella ruminicola et rel.        Prevotella tannerae et rel. 
##                              46918                             744898 
##                 Tannerella et rel.           Uncultured Bacteroidetes 
##                             465584                              44595
```


### Merging operations for phyloseq objects


```r
# merge_phyloseq()
# merge_taxa()
# tax_glom()
```


### Rarification


```r
pseq.rarified <- rarefy_even_depth(pseq)
```

