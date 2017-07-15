<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - diversity}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Global Ecosystem State Variables
--------------------------------

Load example data:

    library(microbiome)
    data(dietswap)
    pseq <- dietswap

### Global indicators

A comprehensive list of global indicators of the ecosystem state can be
obtained as follows. This includes various measures of richness,
evenness, diversity, dominance, and rarity with default parameters. See
the individual functions for more options regarding parameter tuning.

    tab <- global(pseq, index = "all")
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">richness_0</th>
<th align="right">richness_20</th>
<th align="right">richness_50</th>
<th align="right">richness_80</th>
<th align="right">diversities_inverse_simpson</th>
<th align="right">diversities_gini_simpson</th>
<th align="right">diversities_shannon</th>
<th align="right">diversities_fisher</th>
<th align="right">diversities_coverage</th>
<th align="right">evenness_camargo</th>
<th align="right">evenness_pielou</th>
<th align="right">evenness_simpson</th>
<th align="right">evenness_evar</th>
<th align="right">evenness_bulla</th>
<th align="right">dominance_dbp</th>
<th align="right">dominance_dmn</th>
<th align="right">dominance_absolute</th>
<th align="right">dominance_relative</th>
<th align="right">dominance_simpson</th>
<th align="right">dominance_core_abundance</th>
<th align="right">dominance_gini</th>
<th align="right">rarity_log_modulo_skewness</th>
<th align="right">rarity_low_abundance</th>
<th align="right">rarity_noncore_abundance</th>
<th align="right">rarity_rare_abundance</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">112</td>
<td align="right">104</td>
<td align="right">61</td>
<td align="right">27</td>
<td align="right">7.562984</td>
<td align="right">0.8677771</td>
<td align="right">2.942723</td>
<td align="right">12.16148</td>
<td align="right">4</td>
<td align="right">0.1378045</td>
<td align="right">0.6045614</td>
<td align="right">0.0581768</td>
<td align="right">0.0736465</td>
<td align="right">0.2925991</td>
<td align="right">0.3279166</td>
<td align="right">0.4296966</td>
<td align="right">175035</td>
<td align="right">0.3279166</td>
<td align="right">0.1322229</td>
<td align="right">0.9274756</td>
<td align="right">0.8621955</td>
<td align="right">2.059086</td>
<td align="right">0.0291825</td>
<td align="right">0.0150193</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">118</td>
<td align="right">110</td>
<td align="right">68</td>
<td align="right">39</td>
<td align="right">8.105283</td>
<td align="right">0.8766237</td>
<td align="right">2.824184</td>
<td align="right">11.11824</td>
<td align="right">3</td>
<td align="right">0.1159349</td>
<td align="right">0.5802083</td>
<td align="right">0.0623483</td>
<td align="right">0.0722394</td>
<td align="right">0.2506029</td>
<td align="right">0.2428268</td>
<td align="right">0.4655585</td>
<td align="right">323085</td>
<td align="right">0.2428268</td>
<td align="right">0.1233763</td>
<td align="right">0.9328351</td>
<td align="right">0.8840651</td>
<td align="right">2.058747</td>
<td align="right">0.0302304</td>
<td align="right">0.0350443</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">113</td>
<td align="right">104</td>
<td align="right">71</td>
<td align="right">37</td>
<td align="right">4.292701</td>
<td align="right">0.7670464</td>
<td align="right">2.409584</td>
<td align="right">10.80073</td>
<td align="right">2</td>
<td align="right">0.0919433</td>
<td align="right">0.4950318</td>
<td align="right">0.0330208</td>
<td align="right">0.0608332</td>
<td align="right">0.2233591</td>
<td align="right">0.4593873</td>
<td align="right">0.5602856</td>
<td align="right">837328</td>
<td align="right">0.4593873</td>
<td align="right">0.2329536</td>
<td align="right">0.9513098</td>
<td align="right">0.9080567</td>
<td align="right">2.056009</td>
<td align="right">0.0341229</td>
<td align="right">0.0095056</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">114</td>
<td align="right">106</td>
<td align="right">73</td>
<td align="right">30</td>
<td align="right">7.937365</td>
<td align="right">0.8740136</td>
<td align="right">2.994672</td>
<td align="right">11.62450</td>
<td align="right">4</td>
<td align="right">0.1433967</td>
<td align="right">0.6152338</td>
<td align="right">0.0610567</td>
<td align="right">0.0692447</td>
<td align="right">0.2829995</td>
<td align="right">0.3229230</td>
<td align="right">0.3956421</td>
<td align="right">269963</td>
<td align="right">0.3229230</td>
<td align="right">0.1259864</td>
<td align="right">0.8617545</td>
<td align="right">0.8566033</td>
<td align="right">2.054981</td>
<td align="right">0.0349690</td>
<td align="right">0.0370659</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">112</td>
<td align="right">104</td>
<td align="right">65</td>
<td align="right">25</td>
<td align="right">3.172320</td>
<td align="right">0.6847733</td>
<td align="right">2.108225</td>
<td align="right">11.32472</td>
<td align="right">1</td>
<td align="right">0.0790664</td>
<td align="right">0.4331198</td>
<td align="right">0.0244025</td>
<td align="right">0.0709762</td>
<td align="right">0.2066856</td>
<td align="right">0.5448817</td>
<td align="right">0.6315785</td>
<td align="right">596658</td>
<td align="right">0.5448817</td>
<td align="right">0.3152267</td>
<td align="right">0.9533809</td>
<td align="right">0.9209336</td>
<td align="right">2.060330</td>
<td align="right">0.0443278</td>
<td align="right">0.0107751</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">112</td>
<td align="right">105</td>
<td align="right">73</td>
<td align="right">30</td>
<td align="right">2.955214</td>
<td align="right">0.6616150</td>
<td align="right">2.073329</td>
<td align="right">11.18672</td>
<td align="right">1</td>
<td align="right">0.0811396</td>
<td align="right">0.4259505</td>
<td align="right">0.0227324</td>
<td align="right">0.0685637</td>
<td align="right">0.2084814</td>
<td align="right">0.5692406</td>
<td align="right">0.6428424</td>
<td align="right">709407</td>
<td align="right">0.5692406</td>
<td align="right">0.3383850</td>
<td align="right">0.9255292</td>
<td align="right">0.9188604</td>
<td align="right">2.060149</td>
<td align="right">0.0378701</td>
<td align="right">0.0170730</td>
<td align="right">0</td>
</tr>
</tbody>
</table>

### Alpha diversity

This returns a table with selected diversity indicators. See a separate
page on [Beta diversity](Betadiversity.html).

    tab <- diversities(pseq, index = "all")
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">inverse_simpson</th>
<th align="right">gini_simpson</th>
<th align="right">shannon</th>
<th align="right">fisher</th>
<th align="right">coverage</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">7.562984</td>
<td align="right">0.8677771</td>
<td align="right">2.942723</td>
<td align="right">12.16148</td>
<td align="right">4</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">8.105283</td>
<td align="right">0.8766237</td>
<td align="right">2.824184</td>
<td align="right">11.11824</td>
<td align="right">3</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">4.292701</td>
<td align="right">0.7670464</td>
<td align="right">2.409584</td>
<td align="right">10.80073</td>
<td align="right">2</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">7.937365</td>
<td align="right">0.8740136</td>
<td align="right">2.994672</td>
<td align="right">11.62450</td>
<td align="right">4</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">3.172320</td>
<td align="right">0.6847733</td>
<td align="right">2.108225</td>
<td align="right">11.32472</td>
<td align="right">1</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">2.955214</td>
<td align="right">0.6616150</td>
<td align="right">2.073329</td>
<td align="right">11.18672</td>
<td align="right">1</td>
</tr>
</tbody>
</table>

### Richness

This returns observed richness with given detection threshold(s).

    tab <- richness(pseq)
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">0</th>
<th align="right">20</th>
<th align="right">50</th>
<th align="right">80</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">112</td>
<td align="right">104</td>
<td align="right">61</td>
<td align="right">27</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">118</td>
<td align="right">110</td>
<td align="right">68</td>
<td align="right">39</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">113</td>
<td align="right">104</td>
<td align="right">71</td>
<td align="right">37</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">114</td>
<td align="right">106</td>
<td align="right">73</td>
<td align="right">30</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">112</td>
<td align="right">104</td>
<td align="right">65</td>
<td align="right">25</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">112</td>
<td align="right">105</td>
<td align="right">73</td>
<td align="right">30</td>
</tr>
</tbody>
</table>

### Dominance

The dominance index refers to the abundance of the most abundant
species. Various dominance indices are available (see the function help
for a list of options).

    # Absolute abundances for the single most abundant taxa in each sample
    tab <- dominance(pseq, index = "all")
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">dbp</th>
<th align="right">dmn</th>
<th align="right">absolute</th>
<th align="right">relative</th>
<th align="right">simpson</th>
<th align="right">core_abundance</th>
<th align="right">gini</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">0.3279166</td>
<td align="right">0.4296966</td>
<td align="right">175035</td>
<td align="right">0.3279166</td>
<td align="right">0.1322229</td>
<td align="right">0.9274756</td>
<td align="right">0.8621955</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">0.2428268</td>
<td align="right">0.4655585</td>
<td align="right">323085</td>
<td align="right">0.2428268</td>
<td align="right">0.1233763</td>
<td align="right">0.9328351</td>
<td align="right">0.8840651</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">0.4593873</td>
<td align="right">0.5602856</td>
<td align="right">837328</td>
<td align="right">0.4593873</td>
<td align="right">0.2329536</td>
<td align="right">0.9513098</td>
<td align="right">0.9080567</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">0.3229230</td>
<td align="right">0.3956421</td>
<td align="right">269963</td>
<td align="right">0.3229230</td>
<td align="right">0.1259864</td>
<td align="right">0.8617545</td>
<td align="right">0.8566033</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">0.5448817</td>
<td align="right">0.6315785</td>
<td align="right">596658</td>
<td align="right">0.5448817</td>
<td align="right">0.3152267</td>
<td align="right">0.9533809</td>
<td align="right">0.9209336</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">0.5692406</td>
<td align="right">0.6428424</td>
<td align="right">709407</td>
<td align="right">0.5692406</td>
<td align="right">0.3383850</td>
<td align="right">0.9255292</td>
<td align="right">0.9188604</td>
</tr>
</tbody>
</table>

### Rarity and low abundance

The rarity indices quantify the concentration of rare or low abundance
taxa. Various rarity indices are available (see the function help for a
list of options).

    tab <- rarity(pseq, index = "all")
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">log_modulo_skewness</th>
<th align="right">low_abundance</th>
<th align="right">noncore_abundance</th>
<th align="right">rare_abundance</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">2.059086</td>
<td align="right">0.0291825</td>
<td align="right">0.0150193</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">2.058747</td>
<td align="right">0.0302304</td>
<td align="right">0.0350443</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">2.056009</td>
<td align="right">0.0341229</td>
<td align="right">0.0095056</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">2.054981</td>
<td align="right">0.0349690</td>
<td align="right">0.0370659</td>
<td align="right">0</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">2.060330</td>
<td align="right">0.0443278</td>
<td align="right">0.0107751</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">2.060149</td>
<td align="right">0.0378701</td>
<td align="right">0.0170730</td>
<td align="right">0</td>
</tr>
</tbody>
</table>

### Coverage

The coverage index gives the number of groups needed to have a given
proportion of the ecosystem occupied (by default 0.5 ie 50%).

    tab <- coverage(pseq, threshold = 0.5)
    kable(head(tab))

### Core abundance

The core\_abundance function refers to the relative proportion of the
core species. Non-core abundance provides the complement (1-x; see
noncore\_abundance).

    tab <- core_abundance(pseq, detection = .1/100, prevalence = 50/100)
    kable(head(tab))

<table>
<tbody>
<tr class="odd">
<td align="left">Sample-1</td>
<td align="right">0.9665105</td>
</tr>
<tr class="even">
<td align="left">Sample-2</td>
<td align="right">0.9591670</td>
</tr>
<tr class="odd">
<td align="left">Sample-3</td>
<td align="right">0.9720652</td>
</tr>
<tr class="even">
<td align="left">Sample-4</td>
<td align="right">0.9372630</td>
</tr>
<tr class="odd">
<td align="left">Sample-5</td>
<td align="right">0.9851720</td>
</tr>
<tr class="even">
<td align="left">Sample-6</td>
<td align="right">0.9703755</td>
</tr>
</tbody>
</table>

### Gini index

Gini index is a common measure for inequality in economical income. The
inverse gini index (1/x) can also be used as a community diversity
measure.

    tab <- inequality(pseq)
    kable(head(tab))

<table>
<tbody>
<tr class="odd">
<td align="left">Sample-1</td>
<td align="right">0.8621955</td>
</tr>
<tr class="even">
<td align="left">Sample-2</td>
<td align="right">0.8840651</td>
</tr>
<tr class="odd">
<td align="left">Sample-3</td>
<td align="right">0.9080567</td>
</tr>
<tr class="even">
<td align="left">Sample-4</td>
<td align="right">0.8566033</td>
</tr>
<tr class="odd">
<td align="left">Sample-5</td>
<td align="right">0.9209336</td>
</tr>
<tr class="even">
<td align="left">Sample-6</td>
<td align="right">0.9188604</td>
</tr>
</tbody>
</table>

### Evenness

Various evenness measures are also available.

    tab <- evenness(pseq, "all")
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">camargo</th>
<th align="right">pielou</th>
<th align="right">simpson</th>
<th align="right">evar</th>
<th align="right">bulla</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">0.1378045</td>
<td align="right">0.6045614</td>
<td align="right">0.0581768</td>
<td align="right">0.0736465</td>
<td align="right">0.2925991</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">0.1159349</td>
<td align="right">0.5802083</td>
<td align="right">0.0623483</td>
<td align="right">0.0722394</td>
<td align="right">0.2506029</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">0.0919433</td>
<td align="right">0.4950318</td>
<td align="right">0.0330208</td>
<td align="right">0.0608332</td>
<td align="right">0.2233591</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">0.1433967</td>
<td align="right">0.6152338</td>
<td align="right">0.0610567</td>
<td align="right">0.0692447</td>
<td align="right">0.2829995</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">0.0790664</td>
<td align="right">0.4331198</td>
<td align="right">0.0244025</td>
<td align="right">0.0709762</td>
<td align="right">0.2066856</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">0.0811396</td>
<td align="right">0.4259505</td>
<td align="right">0.0227324</td>
<td align="right">0.0685637</td>
<td align="right">0.2084814</td>
</tr>
</tbody>
</table>
