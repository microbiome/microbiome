<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - Preprocessing}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Processing taxonomic profiling data
-----------------------------------

Instructions to manipulate microbiome data sets using tools from the
[phyloseq package](http://joey711.github.io/phyloseq/) and some
extensions from the [microbiome
package](https://github.com/microbiome/microbiome), including
subsetting, aggregating and filtering.

Load example data:

    library(phyloseq)
    library(microbiome)

    data(atlas1006)   # Load the data
    pseq <- core(subset_samples(atlas1006, nationality == "EasternEurope"), detection = 10^2, prevalence = 50/100) # Rename the data and pick subset for faster examples

### Retrieving data elements from a phyloseq object

A phyloseq object contains OTU table (taxa abundances), sample metadata,
taxonomy table (mapping between OTUs and higher-level taxonomic
classifications), and phylogenetic tree (relations between the taxa).
Some of these are optional.

Pick metadata as data.frame:

    meta <- meta(pseq)

Taxonomy table:

    taxonomy <- tax_table(pseq)

Abundances for taxonomic groups ('OTU table') as a TaxaxSamples matrix:

    # Absolute abundances
    otu.absolute <- abundances(pseq)

    # Relative abundances
    otu.relative <- abundances(pseq, "compositional")

Melting phyloseq data for easier plotting:

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
<td>597</td>
<td align="left">Escherichia coli et rel.</td>
<td align="left">Sample-910</td>
<td align="right">179473</td>
<td align="right">56</td>
<td align="left">male</td>
<td align="left">EasternEurope</td>
<td align="left">NA</td>
<td align="left">27</td>
<td align="right">5.51</td>
<td align="left">NA</td>
<td align="left">910</td>
<td align="right">0</td>
<td align="left">Sample-910</td>
<td align="left">Proteobacteria</td>
<td align="left">Escherichia coli et rel.</td>
</tr>
<tr class="even">
<td>1211</td>
<td align="left">Subdoligranulum variable at rel.</td>
<td align="left">Sample-911</td>
<td align="right">162402</td>
<td align="right">45</td>
<td align="left">male</td>
<td align="left">EasternEurope</td>
<td align="left">NA</td>
<td align="left">27</td>
<td align="right">5.62</td>
<td align="left">NA</td>
<td align="left">911</td>
<td align="right">0</td>
<td align="left">Sample-911</td>
<td align="left">Clostridium cluster IV</td>
<td align="left">Subdoligranulum variable at rel.</td>
</tr>
<tr class="odd">
<td>1215</td>
<td align="left">Subdoligranulum variable at rel.</td>
<td align="left">Sample-919</td>
<td align="right">144757</td>
<td align="right">64</td>
<td align="left">male</td>
<td align="left">EasternEurope</td>
<td align="left">NA</td>
<td align="left">28</td>
<td align="right">5.47</td>
<td align="left">NA</td>
<td align="left">919</td>
<td align="right">0</td>
<td align="left">Sample-919</td>
<td align="left">Clostridium cluster IV</td>
<td align="left">Subdoligranulum variable at rel.</td>
</tr>
<tr class="even">
<td>1201</td>
<td align="left">Subdoligranulum variable at rel.</td>
<td align="left">Sample-908</td>
<td align="right">123448</td>
<td align="right">53</td>
<td align="left">male</td>
<td align="left">EasternEurope</td>
<td align="left">NA</td>
<td align="left">27</td>
<td align="right">5.72</td>
<td align="left">NA</td>
<td align="left">908</td>
<td align="right">0</td>
<td align="left">Sample-908</td>
<td align="left">Clostridium cluster IV</td>
<td align="left">Subdoligranulum variable at rel.</td>
</tr>
<tr class="odd">
<td>223</td>
<td align="left">Bifidobacterium</td>
<td align="left">Sample-917</td>
<td align="right">109982</td>
<td align="right">43</td>
<td align="left">male</td>
<td align="left">EasternEurope</td>
<td align="left">NA</td>
<td align="left">28</td>
<td align="right">5.80</td>
<td align="left">NA</td>
<td align="left">917</td>
<td align="right">0</td>
<td align="left">Sample-917</td>
<td align="left">Actinobacteria</td>
<td align="left">Bifidobacterium</td>
</tr>
<tr class="even">
<td>1209</td>
<td align="left">Subdoligranulum variable at rel.</td>
<td align="left">Sample-909</td>
<td align="right">97965</td>
<td align="right">64</td>
<td align="left">female</td>
<td align="left">EasternEurope</td>
<td align="left">NA</td>
<td align="left">27</td>
<td align="right">5.66</td>
<td align="left">NA</td>
<td align="left">909</td>
<td align="right">0</td>
<td align="left">Sample-909</td>
<td align="left">Clostridium cluster IV</td>
<td align="left">Subdoligranulum variable at rel.</td>
</tr>
</tbody>
</table>

### Sample operations

Sample names and variables

    head(sample_names(pseq))

    ## [1] "Sample-312" "Sample-907" "Sample-908" "Sample-909" "Sample-910"
    ## [6] "Sample-911"

Total OTU abundance in each sample

    s <- sample_sums(pseq)

Abundance of a given species in each sample

    head(abundances(pseq)["Akkermansia",])

    ## Sample-312 Sample-907 Sample-908 Sample-909 Sample-910 Sample-911 
    ##       3649       7446       1461       2633       1052       2023

Filtering samples

    f1 <- filterfun_sample(topp(0.1))
    taxa <- genefilter_sample(pseq, f1, A = round(0.5 * nsamples(pseq)))

Select a subset by metadata fields:

    pseq.subset <- subset_samples(pseq, nationality == "AFR")

Select a subset by providing sample names:

    # Check sample names for African Females in this phyloseq object
    s <- rownames(subset(meta(pseq), nationality == "AFR" & sex == "Female"))

    # Pick the phyloseq subset with these sample names
    pseq.subset2 <- prune_samples(s, pseq)

Pick samples at the baseline time points only:

    pseq0 <- baseline(pseq)

### Data transformations

The microbiome package provides a wrapper for standard sample/OTU
transforms. For arbitrary transforms, use the transform\_sample\_counts
function in the phyloseq package. Log10 transform is log(1+x) if the
data contains zeroes. Also "Z", "clr", "hellinger", and "shift" are
available as common transformations. Relative abundances (note that the
input data needs to be in absolute scale, not logarithmic!):

    pseq.compositional <- microbiome::transform(pseq, "compositional")

### Variable operations

Sample variable names

    sample_variables(pseq)

    ##  [1] "age"                   "gender"               
    ##  [3] "nationality"           "DNA_extraction_method"
    ##  [5] "project"               "diversity"            
    ##  [7] "bmi_group"             "subject"              
    ##  [9] "time"                  "sample"

Pick values for a given variable

    head(get_variable(pseq, sample_variables(pseq)[1]))

    ## [1] 36 40 53 64 56 45

Assign new fields to metadata

    # Calculate diversity for samples
    div <- global(pseq, index = "shannon")

    # Assign the estimated diversity to sample metadata
    sample_data(pseq)$diversity <- div

### Taxa operations

Number of taxa

    n <- ntaxa(pseq)

Most abundant taxa

    topx <- top_taxa(pseq, n = 10)

Names

    ranks <- rank_names(pseq)  # Taxonomic levels
    taxa  <- taxa(pseq)        # Taxa names at the analysed level

Subset taxa

    pseq.bac <- subset_taxa(pseq, Phylum == "Bacteroidetes")

Prune (select) taxa:

    # List of Genera in the Bacteroideted Phylum
    taxa <- map_levels(NULL, "Phylum", "Genus", pseq)$Bacteroidetes

    # With given taxon names
    ex2 <- prune_taxa(taxa, pseq)

    # Taxa with positive sum across samples
    ex3 <- prune_taxa(taxa_sums(pseq) > 0, pseq)

Filter by user-specified function values (here variance):

    f <- filter_taxa(pseq, function(x) var(x) > 1e-05, TRUE)

List unique phylum-level groups:

    head(get_taxa_unique(pseq, "Phylum"))

    ## [1] "Verrucomicrobia"          "Proteobacteria"          
    ## [3] "Bacteroidetes"            "Clostridium cluster XIVa"
    ## [5] "Clostridium cluster IV"   "Clostridium cluster XI"

Pick the taxa abundances for a given sample:

    samplename <- sample_names(pseq)[[1]]

    # Pick abundances for a particular taxon
    tax.abundances <- abundances(pseq)[, samplename]

### Merging operations

Aggregate taxa to higher taxonomic levels. This is particularly useful
if the phylogenetic tree is missing. When it is available, see
[merge\_samples, merge\_taxa and
tax\_glom](http://joey711.github.io/phyloseq/merge.html)).

Put the less abundant taxa together in the "Other" category:

    pseq2 <- aggregate_taxa(pseq, "Phylum", top = 5) 

Merging phyloseq objects with the phyloseq package:

    merge_phyloseq(pseqA, pseqB)

### Rarification

    pseq.rarified <- rarefy_even_depth(pseq)

### Taxonomy

Convert between taxonomic levels (here from Genus (Akkermansia) to
Phylum (Verrucomicrobia):

    m <- map_levels("Akkermansia", "Genus", "Phylum", tax_table(pseq))
    print(m)

    ## [1] "Verrucomicrobia"

### Metadata

Visualize frequencies of given factor (sex) levels within the indicated
groups (group):

    res <- plot_frequencies(sample_data(pseq), "bmi_group", "gender")
    print(res$plot)

![](Preprocessing_files/figure-markdown_strict/phylogeny-example3-1.png)

    # Retrieving the actual data values:
    kable(head(res$data), digits = 2)

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
<td align="left">lean</td>
<td align="left">female</td>
<td align="right">1</td>
<td align="right">100.00</td>
</tr>
<tr class="even">
<td align="left">NA</td>
<td align="left">female</td>
<td align="right">6</td>
<td align="right">42.86</td>
</tr>
<tr class="odd">
<td align="left">NA</td>
<td align="left">male</td>
<td align="right">8</td>
<td align="right">57.14</td>
</tr>
</tbody>
</table>
