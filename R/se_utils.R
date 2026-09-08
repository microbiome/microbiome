
# Internal helpers for (Tree)SummarizedExperiment support.
#
# These are deliberately not exported. They exist so that the three
# microbiome accessors (abundances, meta, taxa) can serve both
# phyloseq-class and SummarizedExperiment-derived objects, which in turn
# lets every function that goes through those accessors work with either
# backend without further changes.
#
# SummarizedExperiment is in Suggests, not Imports. That is safe: if the
# user holds a SummarizedExperiment object then the package is by
# definition already loaded, and is() on an undefined class simply
# returns FALSE.


# Is x a SummarizedExperiment (this covers SingleCellExperiment and
# TreeSummarizedExperiment, which both extend it)?
.is_se <- function(x) {
    is(x, "SummarizedExperiment")
}


# Is x one of the container classes we know how to extract a taxa x
# samples matrix from? Used at the class gates in functions that also
# accept plain vectors and matrices.
.is_data_object <- function(x) {
    is.phyloseq(x) || .is_se(x)
}


# Pick a taxa x samples abundance matrix from a SummarizedExperiment.
#
# In SummarizedExperiment features are always rows, so unlike the
# phyloseq path there is no taxa_are_rows() ambiguity to resolve.
#
# assay.type = NULL selects "counts" when present and otherwise the
# first assay, which keeps single-assay objects working without the
# caller having to know the assay name.
.se_assay <- function(x, assay.type=NULL) {

    if (!requireNamespace("SummarizedExperiment", quietly=TRUE)) {
        stop("The SummarizedExperiment package is required to use ",
            "SummarizedExperiment objects with microbiome. ",
            "Please install it.")
    }

    an <- SummarizedExperiment::assayNames(x)
    nassay <- length(SummarizedExperiment::assays(x))

    if (nassay == 0) {
        stop("The SummarizedExperiment object contains no assays.")
    }

    if (is.null(assay.type)) {
        if (!is.null(an) && "counts" %in% an) {
            assay.type <- "counts"
        } else {
            assay.type <- 1L
        }
    } else {
        if (!assay.type %in% an) {
            stop("assay.type '", assay.type, "' not found. Available assays: ",
                paste(an, collapse=", "))
        }
    }

    otu <- SummarizedExperiment::assay(x, assay.type)

    # Sparse and DelayedArray backed assays need materializing, and
    # downstream code assumes a base matrix throughout.
    otu <- as.matrix(otu)

    # assay() does not always carry dimnames; take them from the object.
    rownames(otu) <- rownames(x)
    colnames(otu) <- colnames(x)

    otu

}


# Backend-agnostic counterparts of the phyloseq accessors that are used
# for their side effects or their counts rather than for the abundance
# matrix itself.

.ntaxa <- function(x) {
    .deprecate_phyloseq(x)
    if (.is_se(x)) nrow(x) else ntaxa(x)
}

.nsamples <- function(x) {
    .deprecate_phyloseq(x)
    if (.is_se(x)) ncol(x) else nsamples(x)
}

.sample_names <- function(x) {
    .deprecate_phyloseq(x)
    if (.is_se(x)) colnames(x) else sample_names(x)
}


# Subsetting. The character branch is written as a %in% mask rather than
# as x[keep, ] so that the result keeps the original ordering of x,
# matching what phyloseq::prune_taxa / prune_samples do.

.subset_taxa_obj <- function(x, keep) {

    .deprecate_phyloseq(x)

    if (.is_se(x)) {
        if (is.character(keep)) {
            keep <- rownames(x) %in% keep
        }
        return(x[keep, , drop=FALSE])
    }

    prune_taxa(keep, x)

}

.subset_samples_obj <- function(x, keep) {

    .deprecate_phyloseq(x)

    if (.is_se(x)) {
        if (is.character(keep)) {
            keep <- colnames(x) %in% keep
        }
        return(x[, keep, drop=FALSE])
    }

    prune_samples(keep, x)

}


# Replace the sample metadata. The phyloseq branch keeps the existing
# sample_data<- behaviour.
.set_meta <- function(x, df) {

    .deprecate_phyloseq(x)

    if (.is_se(x)) {
        cd <- methods::as(as.data.frame(df), "DataFrame")
        rownames(cd) <- colnames(x)
        SummarizedExperiment::colData(x) <- cd
        return(x)
    }

    sample_data(x) <- df
    x

}


# Write an abundance matrix back into a data object.
#
# For phyloseq the otu_table is replaced in place, which is the
# established behaviour. For SummarizedExperiment the matrix is stored as
# an additional named assay and the existing assays are left untouched,
# following the mia convention (see transformAssay). The caller therefore
# reads the result back with abundances(x, assay.type = name).
.set_abundances <- function(x, mat, name="counts") {

    .deprecate_phyloseq(x)

    if (.is_se(x)) {

        if (!identical(dim(mat), dim(x))) {
            warning("The transformed matrix has dimensions ",
                paste(dim(mat), collapse="x"), " but the object has ",
                paste(dim(x), collapse="x"),
                ". Returning the matrix instead of the object.")
            return(mat)
        }

        # Drop any stray attributes picked up along the way, for example
        # the scaled:center / scaled:scale that base::scale leaves on the
        # Z transform. Constructing an otu_table discards these on the
        # phyloseq path, so strip them here to keep the two backends
        # returning identical matrices.
        for (a in setdiff(names(attributes(mat)), c("dim", "dimnames"))) {
            attr(mat, a) <- NULL
        }

        rownames(mat) <- rownames(x)
        colnames(mat) <- colnames(x)
        SummarizedExperiment::assay(x, name) <- mat

        return(x)

    }

    if (taxa_are_rows(x)) {
        otu_table(x) <- otu_table(mat, taxa_are_rows=TRUE)
    } else {
        otu_table(x) <- otu_table(t(mat), taxa_are_rows=FALSE)
    }

    x

}
