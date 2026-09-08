
# One-time notice that the phyloseq-based methods are on their way out.
#
# The accessors call each other freely (alpha -> richness -> abundances
# and so on), so emitting this on every call would produce dozens of
# messages per user-level call. It is therefore shown once per session.
# Set options(microbiome.phyloseq_deprecation = FALSE) to silence it, or
# reset the flag below to see it again.

.microbiome_state <- new.env(parent=emptyenv())


.deprecate_phyloseq <- function(x) {

    if (!is.phyloseq(x)) {
        return(invisible(NULL))
    }

    if (!isTRUE(getOption("microbiome.phyloseq_deprecation", TRUE))) {
        return(invisible(NULL))
    }

    if (isTRUE(.microbiome_state$phyloseq_notice_shown)) {
        return(invisible(NULL))
    }

    .microbiome_state$phyloseq_notice_shown <- TRUE

    message("You have inserted a phyloseq object. The phyloseq-based ",
        "methods will be deprecated. We recommend switching to the ",
        "(Tree)SummarizedExperiment framework and mia R/Bioconductor ",
        "package.")

    invisible(NULL)

}


# Show the notice again in the current session. Mainly useful in
# examples, tests and vignettes.
.reset_phyloseq_deprecation <- function() {
    .microbiome_state$phyloseq_notice_shown <- NULL
    invisible(NULL)
}
