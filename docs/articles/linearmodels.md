<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - comparisons}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Linear model examples
---------------------

### Continuous variables

Rapid quantification of continuous associations can be done with the
lm\_phyloseq wrapper function. This uses limma model to generate a table
of P-values and effect sizes. No confounding variables taken into
account in this wrapper. See the [limma
homepage](http://bioinf.wehi.edu.au/limma/) for more info.

    library(limma)
    library(microbiome)
    data("atlas1006")
    # Pick RBB extracted samples (r) and baseline time point
    pseq <- subset_samples(atlas1006, DNA_extraction_method == "r" & time == 0)
    source(system.file("extdata/lm_phyloseq.R", package = "microbiome"))
    tab <- lm_phyloseq(pseq, "age")
    kable(head(tab), digits = 3)

<table>
<thead>
<tr class="header">
<th></th>
<th align="right">logFC</th>
<th align="right">AveExpr</th>
<th align="right">t</th>
<th align="right">P.Value</th>
<th align="right">adj.P.Val</th>
<th align="right">B</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td>Clostridium orbiscindens et rel.</td>
<td align="right">0.005</td>
<td align="right">4.122</td>
<td align="right">4.542</td>
<td align="right">0</td>
<td align="right">0.001</td>
<td align="right">0.696</td>
</tr>
<tr class="even">
<td>Eggerthella lenta et rel.</td>
<td align="right">0.003</td>
<td align="right">2.596</td>
<td align="right">4.451</td>
<td align="right">0</td>
<td align="right">0.001</td>
<td align="right">0.306</td>
</tr>
<tr class="odd">
<td>Butyrivibrio crossotus et rel.</td>
<td align="right">0.004</td>
<td align="right">4.039</td>
<td align="right">4.232</td>
<td align="right">0</td>
<td align="right">0.001</td>
<td align="right">-0.604</td>
</tr>
<tr class="even">
<td>Bifidobacterium</td>
<td align="right">-0.007</td>
<td align="right">4.034</td>
<td align="right">-4.166</td>
<td align="right">0</td>
<td align="right">0.001</td>
<td align="right">-0.872</td>
</tr>
<tr class="odd">
<td>Peptococcus niger et rel.</td>
<td align="right">0.003</td>
<td align="right">2.294</td>
<td align="right">3.872</td>
<td align="right">0</td>
<td align="right">0.003</td>
<td align="right">-2.010</td>
</tr>
<tr class="even">
<td>Eubacterium hallii et rel.</td>
<td align="right">0.004</td>
<td align="right">3.661</td>
<td align="right">3.688</td>
<td align="right">0</td>
<td align="right">0.005</td>
<td align="right">-2.682</td>
</tr>
</tbody>
</table>
