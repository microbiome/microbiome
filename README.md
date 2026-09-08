<!--![Banner](https://github.com/microbiome/microbiome/blob/master/vignettes/figure/composition-example4-1.png)-->
<!--[![Follow](https://img.shields.io/twitter/follow/ropengov.svg?style=social)](https://twitter.com/intent/follow?screen_name=ropengov)-->

microbiome R package
==========

<br>

[![Join the chat at https://gitter.im/microbiome/microbiome](https://badges.gitter.im/microbiome/microbiome.svg)](https://gitter.im/microbiome/microbiome?utm_source=badge&utm_medium=badge&utm_campaign=pr-badge&utm_content=badge)
[![Build Status](https://github.com/microbiome/microbiome/actions/workflows/check-bioc-devel.yml/badge.svg)](https://github.com/microbiome/microbiome/actions/workflows/check-bioc-devel.yml)
[![codecov.io](https://codecov.io/github/microbiome/microbiome/coverage.svg?branch=master)](https://codecov.io/github/microbiome/microbiome?branch=master)
[![PRs Welcome][prs-badge]][prs]
[![Watch on GitHub][github-watch-badge]][github-watch]
[![Star on GitHub][github-star-badge]][github-star]
[![install with bioconda](https://img.shields.io/badge/install%20with-bioconda-brightgreen.svg?style=flat-square)](http://bioconda.github.io/recipes/bioconductor-microbiome/README.html)
<!--[![Follow](https://img.shields.io/twitter/follow/antagomir.svg?style=social)](https://twitter.com/intent/follow?screen_name=antagomir)-->
<!--[![Bioconductor](http://bioconductor.org/shields/build/release/bioc/BiocGenerics.svg)](https://bioconductor.org/packages/devel/bioc/html/microbiome.html)-->


<br>


[prs]: http://makeapullrequest.com
[prs-badge]: https://img.shields.io/badge/PRs-welcome-brightgreen.svg?style=flat-square

[github-watch-badge]: https://img.shields.io/github/watchers/microbiome/microbiome.svg?style=social
[github-watch]: https://github.com/microbiome/microbiome/watchers

[github-star-badge]: https://img.shields.io/github/stars/microbiome/microbiome.svg?style=social
[github-star]: https://github.com/microbiome/microbiome/stargazers
[license-badge]: https://img.shields.io/npm/l/microbiome.svg?style=flat-square
[license]: https://github.com/microbiome/microbiome/blob/master/LICENSE
[microbiome]: https://github.com/microbiome/microbiome



**NOTE** While we continue to maintain this R package, the development
has been discontinued as we have shifted to supporting methods
development based on the new TreeSummarizedExperiment data container,
which provides added capabilities for multi-omics data analysis. Check
the [miaverse project](https://microbiome.github.io/) for details.

**We recommend switching from phyloseq to the TreeSummarizedExperiment
based methods** described in the [Orchestrating Microbiome Analysis
(OMA)](https://bioconductor.org/books/release/OMA/) online book, which
is where method development now takes place.

To ease that transition, most functions in this package now accept
`TreeSummarizedExperiment` (and other `SummarizedExperiment`-derived)
objects in addition to `phyloseq` objects. Two things are worth noting.
Functions that manipulate the taxonomy table (`aggregate_taxa`,
`plot_composition`, `map_levels`, `psmelt2`, the tibble utilities and
the `read_*` family) remain phyloseq only. And for a
`TreeSummarizedExperiment`, `transform()` follows the `mia` convention
of storing its result as a new named assay rather than overwriting the
counts, so the result is read back with
`abundances(x, assay.type = "clr")` rather than `abundances(x)`.


Tools for the exploration and analysis of microbiome profiling data sets.

This R package extends the phyloseq data container, and also supports the (Tree)SummarizedExperiment containers. 
We have discontinued the development and shifted to support methods development based on the (Tree)SummarizedExperiment data containers, see [microbiome.github.io](https://microbiome.github.io/) for more details.

### Installation and use

See the package [tutorial](http://microbiome.github.io/tutorials/).

**Kindly cite** as follows: Sudarshan Shetty and Leo Lahti. Journal of Biosciences 44(5):115, 2019. doi: 10.1007/s12038-019-9930-2
					        
						
### Acknowledgements

Main developer: [Leo Lahti](https://github.com/antagomir/)

Main co-authors: Sudarshan Shetty

[Contributors](https://github.com/microbiome/microbiome/graphs/contributors)

Thanks for [@johanneskoester] and [@nick-youngblut] for contributing [Bioconda installation recipe](https://bioconda.github.io/recipes/bioconductor-microbiome/README.html).

The work has been supported by the following bodies:

  * Academy of Finland (grants 256950, 295741, 307127)
  * [University of Turku](http://www.utu.fi/en/Pages/home.aspx), Department of Mathematics and Statistics
  * [Molecular Ecology group](http://www.mib.wur.nl/UK/), Laboratory of Microbiology, Wageningen University, Netherlands

This work extends the independent [phyloseq](https://github.com/joey711/phyloseq) package and data structures for R-based microbiome analysis. 




