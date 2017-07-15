<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - atlas}
  %\usepackage[utf8]{inputenc}
-->
Intestinal microbiota diversity in 1006 western adults
------------------------------------------------------

The data set from [Lahti et al. Nat. Comm. 5:4344,
2014](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html)
has microbiota profiling of 130 genus-like taxa across 1006 normal
western adults from [Data Dryad](http://doi.org/10.5061/dryad.pk75d).
Load the data in R:

    # Download the required R packages and then the HITChip Atlas data set
    library(microbiome)
    data(atlas1006)

Estimate global ecosystem indicators for this data set:

    tab <- global(atlas1006, index = c("shannon", "invsimpson"))

    library(knitr)
    kable(head(tab))

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">diversities_shannon</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Sample-1</td>
<td align="right">3.189726</td>
</tr>
<tr class="even">
<td>Sample-2</td>
<td align="right">3.396135</td>
</tr>
<tr class="odd">
<td>Sample-3</td>
<td align="right">2.866104</td>
</tr>
<tr class="even">
<td>Sample-4</td>
<td align="right">3.058653</td>
</tr>
<tr class="odd">
<td>Sample-5</td>
<td align="right">3.076850</td>
</tr>
<tr class="even">
<td>Sample-6</td>
<td align="right">2.945709</td>
</tr>
</tbody>
</table>
