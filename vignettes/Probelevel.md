### Probe summarization

Summarize (preprocessed) oligo-level data into phylotype level; examples
with simulated data; see [read.profiling](reading) to use your own data.

    library(microbiome)
    data.directory <- system.file("extdata", package = "microbiome")

    # Read oligo-level data (here: simulated example data)
    oligo.matrix.nolog.simulated <- read.profiling(level = "oligo", 
                                       data.dir = data.directory, log10 = FALSE)

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/oligoprofile.tab

    # Read phylogeny map
    # NOTE: use phylogeny.filtered for species/L1/L2 summarization
    # Load phylogeny.info from output directory
    phylogeny.info <- GetPhylogeny("HITChip", "filtered")

    ## Reading /home/antagomir/R/x86_64-pc-linux-gnu-library/3.1/microbiome/extdata/phylogeny.filtered.tab

    # Summarize oligos into higher level phylotypes
    dat <- summarize.probesets(phylogeny.info, 
                     log10(oligo.matrix.nolog.simulated), 
                     "frpa", "species", verbose = TRUE, 
                      phylotype.rm.list("HITChip"))

    ## Summarizing through species...
    ## Loading pre-calculated preprocessing parameters

### Retrieve probe-level data

Get oligos for each probeset:

    sets <- retrieve.probesets(phylogeny.info, level = "species", name = NULL)

Get probeset data matrix/matrices:

    set <- get.probeset("Actinomyces naeslundii", "species", 
                 phylogeny.info, oligo.matrix.nolog.simulated, log10 = TRUE)
