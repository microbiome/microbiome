<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - Preprocessing}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Preprocessing taxonomic profiling data
--------------------------------------

Here we show how to manipulate microbiome data sets using tools from the
[phyloseq package](http://joey711.github.io/phyloseq/), including
subsetting, aggregating and filtering.

### Picking data from phyloseq

Assuming your data ('pseq' below) is in the phyloseq format, many
standard tools are available.

Sample metadata:

    library(phyloseq)
    library(microbiome)
    data("atlas1006") 
    pseq <- atlas1006
    meta <- sample_data(pseq)

Taxonomy table:

    tax.table <- tax_table(pseq)

Abundances for taxonomic groups ('OTU table') as a taxa x samples
matrix:

    otu <- taxa_abundances(pseq)

Melted data for plotting:

    df <- psmelt(pseq)
    kable(head(df))

<table>
<thead>
<tr class="header">
<th></th>
<th align="left">OTU</th>
<th align="left">Sample</th>
<th align="right">Abundance</th>
<th align="right">age</th>
<th align="left">gender</th>
<th align="left">nationality</th>
<th align="left">DNA_extraction_method</th>
<th align="left">project</th>
<th align="right">diversity</th>
<th align="left">bmi_group</th>
<th align="left">subject</th>
<th align="right">time</th>
<th align="left">sample</th>
<th align="left">Phylum</th>
<th align="left">Genus</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>113110</td>
<td align="left">Prevotella melaninogenica et rel.</td>
<td align="left">Sample-448</td>
<td align="right">944002</td>
<td align="right">54</td>
<td align="left">female</td>
<td align="left">CentralEurope</td>
<td align="left">o</td>
<td align="left">18</td>
<td align="right">5.98</td>
<td align="left">lean</td>
<td align="left">448</td>
<td align="right">0</td>
<td align="left">Sample-448</td>
<td align="left">Bacteroidetes</td>
<td align="left">Prevotella melaninogenica et rel.</td>
</tr>
<tr class="even">
<td>113015</td>
<td align="left">Prevotella melaninogenica et rel.</td>
<td align="left">Sample-360</td>
<td align="right">902034</td>
<td align="right">45</td>
<td align="left">female</td>
<td align="left">CentralEurope</td>
<td align="left">o</td>
<td align="left">13</td>
<td align="right">5.49</td>
<td align="left">severeobese</td>
<td align="left">360</td>
<td align="right">0</td>
<td align="left">Sample-360</td>
<td align="left">Bacteroidetes</td>
<td align="left">Prevotella melaninogenica et rel.</td>
</tr>
<tr class="odd">
<td>112747</td>
<td align="left">Prevotella melaninogenica et rel.</td>
<td align="left">Sample-190</td>
<td align="right">862870</td>
<td align="right">34</td>
<td align="left">female</td>
<td align="left">CentralEurope</td>
<td align="left">r</td>
<td align="left">7</td>
<td align="right">6.06</td>
<td align="left">lean</td>
<td align="left">190</td>
<td align="right">0</td>
<td align="left">Sample-190</td>
<td align="left">Bacteroidetes</td>
<td align="left">Prevotella melaninogenica et rel.</td>
</tr>
<tr class="even">
<td>113109</td>
<td align="left">Prevotella melaninogenica et rel.</td>
<td align="left">Sample-743</td>
<td align="right">852350</td>
<td align="right">52</td>
<td align="left">male</td>
<td align="left">US</td>
<td align="left">NA</td>
<td align="left">19</td>
<td align="right">5.21</td>
<td align="left">obese</td>
<td align="left">743</td>
<td align="right">0</td>
<td align="left">Sample-743</td>
<td align="left">Bacteroidetes</td>
<td align="left">Prevotella melaninogenica et rel.</td>
</tr>
<tr class="odd">
<td>112944</td>
<td align="left">Prevotella melaninogenica et rel.</td>
<td align="left">Sample-366</td>
<td align="right">851147</td>
<td align="right">52</td>
<td align="left">female</td>
<td align="left">CentralEurope</td>
<td align="left">o</td>
<td align="left">15</td>
<td align="right">5.63</td>
<td align="left">obese</td>
<td align="left">366</td>
<td align="right">0</td>
<td align="left">Sample-366</td>
<td align="left">Bacteroidetes</td>
<td align="left">Prevotella melaninogenica et rel.</td>
</tr>
<tr class="even">
<td>113639</td>
<td align="left">Prevotella melaninogenica et rel.</td>
<td align="left">Sample-375</td>
<td align="right">844482</td>
<td align="right">45</td>
<td align="left">female</td>
<td align="left">CentralEurope</td>
<td align="left">o</td>
<td align="left">16</td>
<td align="right">5.64</td>
<td align="left">severeobese</td>
<td align="left">375</td>
<td align="right">0</td>
<td align="left">Sample-375</td>
<td align="left">Bacteroidetes</td>
<td align="left">Prevotella melaninogenica et rel.</td>
</tr>
</tbody>
</table>

### Standard data processing operations

Let us look at [a two-week diet swap
study](http://dx.doi.org/10.1038/ncomms7342) between western (USA) and
traditional (rural Africa) diets including microbiota profiling:

    library(microbiome)
    data("dietswap")
    pseq <- dietswap

### Sample operations

Sample names and variables

    head(sample_names(pseq))

    ## [1] "Sample-1" "Sample-2" "Sample-3" "Sample-4" "Sample-5" "Sample-6"

Sample sums

    head(sample_sums(pseq))

    ## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
    ##   533779  1330516  1822706   835998  1095023  1246234

Abundance of a given species in each sample

    head(get_sample(pseq, taxa_names(pseq)[1]))

    ## Sample-1 Sample-2 Sample-3 Sample-4 Sample-5 Sample-6 
    ##       11       67       21       42       16       20

Filter samples

    f1 <- filterfun_sample(topp(0.1))
    taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))

Select samples by specific metadata fields:

    pseq.subset <- subset_samples(pseq, nationality == "AFR")

### Data transformations

The microbiome package provides a wrapper for standard sample/OTU
transformations. For arbitrary transformations, use the
transform\_sample\_counts function in the phyloseq package.

Log10 transformation (log(1+x) if the data contains zeroes)

    pseq.log <- transform_phyloseq(pseq, "log10")

Z transformation:

    pseq.zotu <- transform_phyloseq(pseq, "Z", "OTU")

Relative abundances (the input data needs to be in absolute scale, not
logarithmic!):

    pseq1 <- transform_phyloseq(pseq, "compositional", "OTU")
    pseq2 <- transform_sample_counts(pseq, function(x) x/sum(x))

### Variable operations

Sample variable names

    sample_variables(pseq)

    ## [1] "subject"                "sex"                   
    ## [3] "nationality"            "group"                 
    ## [5] "sample"                 "timepoint"             
    ## [7] "timepoint.within.group" "bmi_group"

Pick variable values for a given variable

    head(get_variable(pseq, sample_variables(pseq)[1]))

    ## [1] byn nms olt pku qjy riv
    ## 38 Levels: azh azl byn byu cxj dwc dwk eve fua fud gtd gty hsf irh ... zaq

    # .. or assigning fields to metadata:
    # sample_data(GP)$human <- ..

### Taxa operations

Number of taxa

    n <- ntaxa(pseq)

Names

    ranks <- rank_names(pseq)
    taxa <- taxa_names(pseq)

Prune taxa:

    taxa <- map_levels(NULL, "Phylum", "Genus", pseq)$Bacteroidetes

    # With given taxon names
    ex2 <- prune_taxa(taxa, pseq)

    # Taxa with positive sum across samples
    ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)

Subset taxa:

    pseq <- subset_taxa(pseq, Phylum == "Bacteroidetes")

Filter by user-specified function values (here variance):

    f <- filter_taxa(pseq, function(x) var(x) > 1e-05, TRUE)

List unique phylum-level groups:

    head(get_taxa_unique(pseq, "Phylum"))

    ## [1] "Bacteroidetes"

Pick detected taxa by sample name:

    samplename <- sample_names(pseq)[[1]]
    tax.abundances <- get_taxa(pseq, samplename)

Taxa sums

    head(taxa_sums(pseq))

    ##               Allistipes et rel.     Bacteroides fragilis et rel. 
    ##                          3513027                          2539567 
    ## Bacteroides intestinalis et rel.       Bacteroides ovatus et rel. 
    ##                           199684                          1516522 
    ##     Bacteroides plebeius et rel.  Bacteroides splachnicus et rel. 
    ##                           596972                           833871

### Merging operations

Aggregate taxa to higher taxonomic levels. This is particularly useful
if the phylogenetic tree is missing. When it is available, see
[merge\_samples, merge\_taxa and
tax\_glom](http://joey711.github.io/phyloseq/merge.html))

    pseq2 <- summarize_taxa(pseq, "Phylum") 

Merging phyloseq objects

    merge_phyloseq(pseqA, pseqB)

### Rarification

    pseq.rarified <- rarefy_even_depth(pseq)

### Taxonomy

Convert between taxonomic levels (here from Genus (Akkermansia) to
Phylum (Verrucomicrobia)):

    data(atlas1006)
    pseq <- atlas1006
    m <- map_levels("Akkermansia", "Genus", "Phylum", tax_table(pseq))
    print(m)

    ## [1] "Verrucomicrobia"

### Metadata

Visualize frequencies of given factor (sex) levels within the indicated
groups (group):

    data(dietswap)
    res <- plot_frequencies(sample_data(dietswap), "group", "sex")
    print(res$plot)

![](Preprocessing_files/figure-markdown_strict/phylogeny-example3-1.png)

    kable(res$data, digits = 2)

<table>
<thead>
<tr class="header">
<th align="left">Groups</th>
<th align="left">Factor</th>
<th align="right">n</th>
<th align="right">pct</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">DI</td>
<td align="left">Female</td>
<td align="right">34</td>
<td align="right">47.22</td>
</tr>
<tr class="even">
<td align="left">DI</td>
<td align="left">Male</td>
<td align="right">38</td>
<td align="right">52.78</td>
</tr>
<tr class="odd">
<td align="left">ED</td>
<td align="left">Female</td>
<td align="right">34</td>
<td align="right">45.33</td>
</tr>
<tr class="even">
<td align="left">ED</td>
<td align="left">Male</td>
<td align="right">41</td>
<td align="right">54.67</td>
</tr>
<tr class="odd">
<td align="left">HE</td>
<td align="left">Female</td>
<td align="right">34</td>
<td align="right">45.33</td>
</tr>
<tr class="even">
<td align="left">HE</td>
<td align="left">Male</td>
<td align="right">41</td>
<td align="right">54.67</td>
</tr>
</tbody>
</table>
