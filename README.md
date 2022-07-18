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


Tools for the exploration and analysis of microbiome profiling data
sets.

This R package extends the phyloseq data container. The package is actively maintened but we have discontinued the development and shifted to support methods development based on the (Tree)SummarizedExperiment data containers, see [microbiome.github.io](https://microbiome.github.io/) for more details.

### Installation and use

See the package [tutorial](http://microbiome.github.io/tutorials/).

**Kindly cite** as follows: "Leo Lahti, Sudarshan Shetty [et al.](https://github.com/microbiome/microbiome/graphs/contributors) ([Bioconductor, 2017](https://bioconductor.org/packages/devel/bioc/html/microbiome.html)). Tools for microbiome analysis in R. Microbiome package version 1.17.43. URL: [http://microbiome.github.com/microbiome](http://microbiome.github.com/microbiome). See also the relevant references listed in the manual page of each function. 

### Contribute

Contributions and feedback are very welcome:

  * [Issue Tracker](https://github.com/microbiome/microbiome/issues) 
  * [Pull requests](https://github.com/microbiome/microbiome/)
  * Subscribe to the [mailing list](https://groups.google.com/forum/#!forum/microbiome-devel) (microbiome-devel@googlegroups.com)
  * [Gitter chat room](https://gitter.im/microbiome)
  * [Star us on the Github page](https://github.com/microbiome/microbiome)


### Publications using the microbiome package

Below some publications that utilize the tools implemented in this package. The list of publications is not exhaustive. Let us know if you know of further publications using the microbiome package; we are collecting these on the website.

[Intestinal microbiome landscaping: Insight in community assemblage and implications for microbial modulation strategies](https://academic.oup.com/femsre/article/doi/10.1093/femsre/fuw045/2979411/Intestinal-microbiome-landscaping-insight-in#58802539). Shetty S, Hugenholtz F, Lahti L, Smidt H, de Vos WM, Danchin A. _FEMS Microbiology Reviews_ fuw045, 2017.

[Metagenomics meets time series analysis: unraveling microbial community dynamics](http://dx.doi.org/10.1016/j.mib.2015.04.004) Faust K, Lahti L, Gonze D, de Vos WM, Raes J. _Current Opinion in Microbiology_ 15:56-66 2015.

[Tipping elements in the human intestinal ecosystem](http://www.nature.com/ncomms/2014/140708/ncomms5344/full/ncomms5344.html) Lahti L, Salojärvi J, Salonen A, Scheffer M, de Vos WM. _Nature Communications_ 5:4344, 2014. 

[Fat, Fiber and Cancer Risk in African, Americans and Rural Africans](http://www.nature.com/ncomms/2015/150428/ncomms7342/full/ncomms7342.html)  O’Keefe S, Li JV, Lahti L, Ou J, Carbonero F, Mohammed K, Posma JM, Kinross J, Wahl E, Ruder E, Vipperla K, Naidoo V, Mtshali L, Tims S, Puylaert PGB, DeLany J, Krasinskas A, Benefiel AC, Kaseb HO, Newton K, Nicholson JK, de Vos WM, Gaskins HR, Zoetendal EG. _Nature Communications_ 6:6342, 2015.

[Associations between the human intestinal microbiota, Lactobacillus rhamnosus GG and serum lipids indicated by integrated analysis of high-throughput profiling data](http://dx.doi.org/10.7717/peerj.32) Lahti L, Salonen A, Kekkonen RA, Salojärvi J, Jalanka-Tuovinen J, Palva A, Orešič M, de Vos WM. _PeerJ_ 1:e32, 2013.

[The adult intestinal core microbiota is determined by analysis depth and health status](http://onlinelibrary.wiley.com/doi/10.1111/j.1469-0691.2012.03855.x/abstract) Salonen A, Salojärvi J, Lahti L, and de Vos WM. _Clinical Microbiology and Infection_ 18(S4):16 20, 2012. 


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




