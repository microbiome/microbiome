context("TreeSummarizedExperiment support")

# These tests check that the functions which go through the microbiome
# accessors (abundances, meta, taxa) return identical results for a
# phyloseq object and for the equivalent TreeSummarizedExperiment.

skip_if_no_tse <- function() {
    skip_if_not_installed("mia")
    skip_if_not_installed("TreeSummarizedExperiment")
}

test_that("accessors agree between phyloseq and TreeSE", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    expect_equal(abundances(dietswap), abundances(tse))
    expect_equal(taxa(dietswap), taxa(tse))
    expect_equal(dim(meta(dietswap)), dim(meta(tse)))
    expect_equal(rownames(meta(dietswap)), rownames(meta(tse)))
})


test_that("assay.type selects the requested assay", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    # Add a second assay and check it is reachable and distinct
    SummarizedExperiment::assay(tse, "doubled") <-
        SummarizedExperiment::assay(tse, "counts") * 2

    expect_equal(abundances(tse, assay.type="counts") * 2,
        abundances(tse, assay.type="doubled"))

    # Default still picks counts when present
    expect_equal(abundances(tse), abundances(tse, assay.type="counts"))

    expect_error(abundances(tse, assay.type="nonexistent"), "not found")
})


test_that("tier A functions agree between phyloseq and TreeSE", {

    skip_if_no_tse()

    data(dietswap)
    ps <- dietswap
    tse <- mia::convertFromPhyloseq(ps)

    # Diversity and dominance
    expect_equal(alpha(ps, index="shannon"), alpha(tse, index="shannon"))
    expect_equal(microbiome::diversity(ps), microbiome::diversity(tse))
    expect_equal(richness(ps), richness(tse))
    expect_equal(evenness(ps), evenness(tse))
    expect_equal(dominance(ps), dominance(tse))
    expect_equal(dominant(ps), dominant(tse))
    expect_equal(rarity(ps), rarity(tse))
    expect_equal(inequality(ps), inequality(tse))

    # Divergence against a common reference
    ref <- apply(abundances(ps), 1, median)
    expect_equal(divergence(ps, ref), divergence(tse, ref))

    # Summaries
    expect_equal(microbiome::coverage(ps), microbiome::coverage(tse))
    expect_equal(readcount(ps), readcount(tse))
    expect_equal(low_abundance(ps), low_abundance(tse))
    expect_equal(log_modulo_skewness(ps), log_modulo_skewness(tse))
    expect_equal(is_compositional(ps), is_compositional(tse))
    expect_equal(overlap(ps), overlap(tse))
    expect_equal(bimodality(ps, method="Sarle.finite.sample"),
        bimodality(tse, method="Sarle.finite.sample"))

    # Core / rare
    expect_equal(prevalence(ps), prevalence(tse))
    expect_equal(core_members(ps, detection=1, prevalence=0.5),
        core_members(tse, detection=1, prevalence=0.5))
    expect_equal(core_abundance(ps), core_abundance(tse))
    expect_equal(rare_members(ps), rare_members(tse))
    expect_equal(rare_abundance(ps), rare_abundance(tse))

    expect_equal(top_taxa(ps, n=10), top_taxa(tse, n=10))
})


test_that("transform stores the result as a new named assay", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    for (method in c("compositional", "clr", "log10", "Z")) {

        res <- microbiome::transform(tse, method)

        expect_true(is(res, "TreeSummarizedExperiment"))

        # The transformed values match the phyloseq result
        expect_equal(
            abundances(res, assay.type=method),
            abundances(microbiome::transform(dietswap, method)))

        # The original counts are left untouched (option a)
        expect_true("counts" %in% SummarizedExperiment::assayNames(res))
        expect_equal(abundances(res), abundances(tse))
    }

    # The assay name is configurable
    res <- microbiome::transform(tse, "compositional", name="relabundance")
    expect_true("relabundance" %in% SummarizedExperiment::assayNames(res))
    expect_equal(abundances(res, assay.type="relabundance"),
        abundances(microbiome::transform(dietswap, "compositional")))
})


test_that("subsetting functions agree between backends", {

    skip_if_no_tse()

    data(dietswap)
    ps <- dietswap
    tse <- mia::convertFromPhyloseq(ps)

    # core / rare filter taxa
    expect_equal(taxa(core(ps, detection=1, prevalence=0.5)),
        taxa(core(tse, detection=1, prevalence=0.5)))
    expect_equal(taxa(rare(ps, detection=1, prevalence=0.5)),
        taxa(rare(tse, detection=1, prevalence=0.5)))

    # remove_taxa / remove_samples
    drop_t <- taxa(ps)[seq_len(5)]
    expect_equal(taxa(remove_taxa(drop_t, ps)), taxa(remove_taxa(drop_t, tse)))

    drop_s <- .sample_names(ps)[seq_len(5)]
    expect_equal(colnames(abundances(remove_samples(drop_s, ps))),
        colnames(abundances(remove_samples(drop_s, tse))))

    # The subset keeps its class
    expect_true(is(core(tse, detection=1, prevalence=0.5),
        "TreeSummarizedExperiment"))
    expect_true(is(remove_samples(drop_s, tse), "TreeSummarizedExperiment"))
})


test_that("time series helpers agree between backends", {

    skip_if_no_tse()

    data(dietswap)
    ps <- dietswap
    # dietswap carries 'timepoint' rather than the 'time' field these
    # functions expect
    sample_data(ps)$time <- sample_data(ps)$timepoint
    tse <- mia::convertFromPhyloseq(ps)

    expect_equal(colnames(abundances(baseline(ps))),
        colnames(abundances(baseline(tse))))
    expect_equal(timesplit(ps), timesplit(tse))
    expect_equal(time_sort(ps), time_sort(tse))
    expect_equal(meta(time_normalize(ps))$time, meta(time_normalize(tse))$time)

    # time_normalize writes metadata back into the object
    expect_true(is(time_normalize(tse), "TreeSummarizedExperiment"))

    set.seed(1); a <- .sample_names(collapse_replicates(ps,
        replicate_fields="subject"))
    set.seed(1); b <- .sample_names(collapse_replicates(tse,
        replicate_fields="subject"))
    expect_equal(sort(a), sort(b))
})


test_that("summarize_phyloseq agrees between backends", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    expect_equal(suppressMessages(summarize_phyloseq(dietswap)),
        suppressMessages(summarize_phyloseq(tse)))
})


test_that("plotting functions accept TreeSE", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    expect_s3_class(
        plot_core(tse, prevalences=seq(.1, 1, .2), detections=3), "ggplot")
    expect_s3_class(
        plot_tipping(tse, taxon=taxa(tse)[[1]], tipping.point=0.5), "ggplot")
    expect_s3_class(
        boxplot_abundance(tse, x="nationality", y=taxa(tse)[[1]]), "ggplot")
    expect_s3_class(
        suppressMessages(boxplot_alpha(tse, x_var="nationality",
            index="shannon")), "ggplot")
    expect_s3_class(plot_atlas(tse, x="nationality", y=taxa(tse)[[1]]),
        "ggplot")
    expect_s3_class(plot_density(tse, variable=taxa(tse)[[1]]), "ggplot")
})


test_that("compositional transform of abundances agrees", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    expect_equal(abundances(dietswap, transform="compositional"),
        abundances(tse, transform="compositional"))
    expect_equal(abundances(dietswap, transform="clr"),
        abundances(tse, transform="clr"))
})


test_that("phyloseq input emits a one-time deprecation notice", {

    data(dietswap)

    # Fires on first phyloseq contact
    microbiome:::.reset_phyloseq_deprecation()
    expect_message(abundances(dietswap), "phyloseq-based methods will be")

    # ... but only once per session, so that the internal accessor calls
    # do not bury the user in repeats
    expect_no_message(abundances(dietswap))
    expect_no_message(meta(dietswap))
    expect_no_message(taxa(dietswap))

    # Honours the option
    microbiome:::.reset_phyloseq_deprecation()
    old <- options(microbiome.phyloseq_deprecation=FALSE)
    on.exit(options(old), add=TRUE)
    expect_no_message(abundances(dietswap))
})


test_that("TreeSE input never emits the phyloseq notice", {

    skip_if_no_tse()

    data(dietswap)
    tse <- mia::convertFromPhyloseq(dietswap)

    microbiome:::.reset_phyloseq_deprecation()
    expect_no_message(abundances(tse))
    expect_no_message(meta(tse))
    expect_no_message(taxa(tse))
    expect_no_message(core(tse, detection=1, prevalence=0.5))
})
