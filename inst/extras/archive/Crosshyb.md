<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - crosshyb}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Visualizing cross-hybridization
-------------------------------

To visualize cross-hybridization between selected taxa on HITChip (or
other chips), use the following scripts. By default the groups with no
cross-hyb are filtered out for clarity. Rows and columns are ordered by
hierarchical clustering. The cross-hyb is shown in percentages, rounded
as indicated by the rounding argument. The percentage indicates which
fraction of probes for a row taxon overlaps with probes of a column
taxon. This is not symmetric if the row and col taxon have a different
total number of probes. For details, see help(plot\_crosshyb).

    library(microbiome, quietly = TRUE)
    library(dplyr)

    # Pick the phylogeny which was used to summarize probes to species level
    tax.table <- get_hitchip_taxonomy("HITChip", "full")

    # Check cross-hyb between all L2 groups
    res <- plot_crosshyb(tax.level = "L2", rounding = 1, show.plot = FALSE, tax.table = tax.table)
        
    # Pick the crosshyb table and figure
    crosshyb.table <- res$data
    p <- res$plot

    # Plot the figure    
    print(p)

![](Crosshyb_files/figure-markdown_strict/chyb-1.png)

    # Organize the Crosshyb table
    suppressMessages(library(dplyr))
    s <- filter(res$data, crosshyb > 0)
    s <- s[rev(order(s$crosshyb)),]
    head(s)

    ##                          Taxon1                         Taxon2  crosshyb
    ## 886    Uncultured Bacteroidetes   Bacteroides plebeius et rel. 100.00000
    ## 244    Uncultured Bacteroidetes             Allistipes et rel.  96.42857
    ## 6539   Uncultured Bacteroidetes             Tannerella et rel.  92.85714
    ## 7043                Leminorella               Yersinia et rel.  91.30435
    ## 780  Bacteroides ovatus et rel.   Bacteroides fragilis et rel.  90.90909
    ## 4703               Burkholderia Oxalobacter formigenes et rel.  90.00000

### Further examples

Investigate species-species cross-hybridization within the Dialister L2
group

    # Select species belonging to Dialister L2 group
    mytaxa <- map_levels("Dialister", from = "L2", to = "species", tax.table)[[1]]

    # Check cross-hyb between Dialister species
    res <- plot_crosshyb(tax.level = "species", selected.taxa = mytaxa, rounding = 0, tax.table = tax.table)

![](Crosshyb_files/figure-markdown_strict/chyb2-1.png)

    # Check the cross-hyb data as well
    library(knitr)
    kable(head(res$data))

<table>
<thead>
<tr class="header">
<th align="left">Taxon1</th>
<th align="left">Taxon2</th>
<th align="right">crosshyb</th>
</tr>
</thead>
<tbody>
<tr class="odd">
<td align="left">Dialister invisus</td>
<td align="left">Dialister invisus</td>
<td align="right">0</td>
</tr>
<tr class="even">
<td align="left">Dialister pneumosintes</td>
<td align="left">Dialister invisus</td>
<td align="right">50</td>
</tr>
<tr class="odd">
<td align="left">Uncultured bacterium clone Eldhufec089</td>
<td align="left">Dialister invisus</td>
<td align="right">100</td>
</tr>
<tr class="even">
<td align="left">Uncultured bacterium clone Eldhufec093</td>
<td align="left">Dialister invisus</td>
<td align="right">100</td>
</tr>
<tr class="odd">
<td align="left">Uncultured bacterium clone Eldhufec096</td>
<td align="left">Dialister invisus</td>
<td align="right">100</td>
</tr>
<tr class="even">
<td align="left">uncultured bacterium MG10</td>
<td align="left">Dialister invisus</td>
<td align="right">40</td>
</tr>
</tbody>
</table>
