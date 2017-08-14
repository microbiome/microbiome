<!--
  %\VignetteEngine{knitr::rmarkdown}
  %\VignetteIndexEntry{microbiome tutorial - potential}
  %\usepackage[utf8]{inputenc}
  %\VignetteEncoding{UTF-8}  
-->
Potential analysis
------------------

Potential analysis (following [Hirota et al. Science, 334,
232-235.](http://www.sciencemag.org/content/334/6053/232.long)) provides
tools to assess how states of an indicator variable vary with respect to
a given background variable.

Load example data:

    library(microbiome)
    data(atlas1006) # From http://doi.org/10.5061/dryad.pk75d
    pseq <- atlas1006 

Assess the relationship between age and microbiome diversity:

    # Pick diversity and age
    diversity <- exp(microbiome::diversity(pseq)$Shannon)
    age <- meta(pseq)$age

    # Run potential analysis
    library(earlywarnings)
    res <- movpotential_ews(diversity, age)

Visualize

    p <- plot_potential(res$res) + xlab("Age") + ylab("Diversity")
    print(p)
