#' @title Bimodality Analysis
#' @description Estimate bimodality scores.
#' @param x A vector, matrix, or a phyloseq object
#' @param method bimodality quantification method
#' ('potential_analysis', 'Sarle.finite.sample', or 'Sarle.asymptotic').
#' If method='all', then a data.frame with all scores is returned.
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iter Bootstrap iterations
#' @param min.density minimum accepted density for a maximum;
#' as a multiple of kernel height
#' @param verbose Verbose
#' @inheritParams potential_analysis
#' @return A list with following elements:
#' \itemize{
#' \item{score}{Fraction of bootstrap samples where multiple modes are
#' observed}
#' \item{nmodes}{The most frequently observed number of modes in
#' bootstrap sampling results.}
#' \item{results}{Full results of potential_analysis for each row of the
#' input matrix.}
#' }
#' @details
#' \itemize{
#' \item{Sarle.finite.sample}{ Coefficient of bimodality for
#' finite sample. See SAS 2012.}
#' \item{Sarle.asymptotic}{ Coefficient of bimodality, used and described
#' in Shade et al. (2014) and Ellison AM (1987).}
#' \item{potential_analysis}{ Repeats potential analysis
#' (Livina et al. 2010) multiple times with bootstrap sampling for
#' each row of the input data (as in Lahti et al. 2014) and returns
#' the bootstrap score.}
#' }
#' 
#' The coefficient lies in (0, 1).
#' 
#' The 'Sarle.asymptotic' version is defined as
#' \deqn{b=(g^2 + 1) / k}.
#' This is coefficient of bimodality from Ellison
#' AM Am. J. Bot. 1987, for microbiome analysis it has been used for
#' instance in Shade et al. 2014.
#' The formula for 'Sarle.finite.sample' (SAS 2012):
#' \deqn{b=\frac{g^2 + 1}{k + (3(n-1)^2)/((n-2)(n-3))}}
#' where n is sample size and 
#' In both formulas, \eqn{g} is sample skewness and \eqn{k} is the kth
#' standardized moment (also called the sample kurtosis, or
#' excess kurtosis).
#'
#' @seealso A classical test of multimodality is provided by \code{dip.test}
#' in the \pkg{DIP} package.
#'
#' @references
#' \itemize{
#' \item{}{Livina et al. (2010). Potential analysis 
#' reveals changing number of climate states during the last 60
#' kyr. \emph{Climate of the Past}, 6, 77-82.}
#' \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.}
#' \item{}{Shade et al. mBio 5(4):e01371-14, 2014.}
#' \item{}{AM Ellison, Am. J. Bot 74:1280-8, 1987.}
#' \item{}{SAS Institute Inc. (2012). SAS/STAT 12.1 user's guide. Cary, NC.}
#' \item{}{To cite the microbiome R package, see citation('microbiome')}
#' }
#' @export
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#' # In practice, use more bootstrap iterations   
#' b <- bimodality(c(rnorm(100, mean=0), rnorm(100, mean=5)),
#'     method = "Sarle.finite.sample", bs.iter=5)
#' # The classical DIP test:
#' # quantifies unimodality. Values range between 0 to 1. 
#' # dip.test(x, simulate.p.value=TRUE, B=200)$statistic
#' # Values less than 0.05 indicate significant deviation from unimodality.
#' # Therefore, to obtain an increasing multimodality score, use
#' # library(diptest)
#' # multimodality.dip <- apply(abundances(pseq), 1,
#' # function (x) {1 - unname(dip.test(x)$p.value)})
#'
#' @keywords utilities
bimodality <- function(x, method="potential_analysis", peak.threshold=1,
    bw.adjust=1, bs.iter=100, min.density=1, verbose=TRUE) {
    
    accepted <- intersect(method, c("potential_analysis",
        "Sarle.finite.sample", "Sarle.asymptotic"))
    
    if (length(method) > 1 || method == "all") {
        method <- accepted
        tab <- NULL
        for (meth in method) {
            b <- bimodality(x, method=meth, peak.threshold, bw.adjust,
                bs.iter, min.density, verbose)
            tab <- cbind(tab, b)
        }
        colnames(tab) <- method
        tab <- as.data.frame(tab)
        return(tab)
    }
    
    
    if (is.vector(x)) {
        
        if (method %in% c("Sarle.finite.sample", "Sarle.asymptotic")) {
            
            s <- bimodality_sarle(x, type=method)
            
        } else if (method == "potential_analysis") {

            if (length(unique(x)) == 1) {                
                s <- 0               
            } else {
                # Shift the data. This does not affect mode detection but
                # avoids errors with nonnegatives.
                s <- multimodality(x, peak.threshold, bw.adjust,
                                    bs.iter, min.density, verbose)$score
            }
        }
        
    } else if (is.matrix(x)) {
        
        s <- apply(x, 1, function(xi) {
            bimodality(xi, method=method, peak.threshold=peak.threshold,
            bw.adjust=bw.adjust, 
                bs.iter=bs.iter, min.density=min.density, verbose=verbose)
        })
        
    } else if (is.phyloseq(x)) {
        
        # Pick the data from phyloseq object
        x <- abundances(x)
        s <- bimodality(x, method=method, peak.threshold=peak.threshold,
        bw.adjust=bw.adjust, 
            bs.iter=bs.iter, min.density=min.density, verbose=verbose)
        
    }
    
    s
    
}





#' @title Multimodality Score
#' @description Multimodality score based on bootstrapped potential analysis.
#' @param x A vector, or data matrix (variables x samples)
#' @param bw.adjust Bandwidth adjustment
#' @param bs.iter Bootstrap iterations
#' @param min.density minimum accepted density for a maximum;
#' as a multiple of kernel height
#' @param verbose Verbose
#' @inheritParams potential_analysis
#' @return A list with following elements: 
#' \itemize{
#' \item{score}{Fraction of bootstrap samples with multiple
#' observed modes}
#' \item{nmodes}{The most frequently observed number of modes
#' in bootstrap}
#' \item{results}{Full results of potential_analysis for each
#' row of the input matrix.}
#' }
#' @details Repeats potential analysis (Livina et al. 2010) multiple times
#' with bootstrap sampling for each row of the input data
#' (as in Lahti et al. 2014) and returns the specified results.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples
#' #data(peerj32)
#' #s <- multimodality(t(peerj32$microbes[, c('Akkermansia', 'Dialister')]))
#' @references
#' \itemize{
#' \item{}{Livina et al. (2010). Potential analysis reveals changing
#' number of climate states during the last 60 kyr.
#' \emph{Climate of the Past}, 6, 77-82.}
#' \item{}{Lahti et al. (2014). Tipping elements of the human intestinal
#' ecosystem. \emph{Nature Communications} 5:4344.}
#' }
#' @keywords utilities
multimodality <- function(x, peak.threshold=1, bw.adjust=1,
    bs.iter=100, min.density=1, verbose=TRUE) {
    
    if (is.vector(x)) {

        # Add small noise to enable robust density estimation
        # (identical values may cause failure)
        x <- x + rnorm(length(x), sd=sd(x)/100)
        m <- potential_analysis(x, bw.adjust=bw.adjust, bs.iter=bs.iter,
                min.density=min.density)
        ret <- list(score=1 - m$unimodality.support,
        modes=m$modes, results=m)
        return(ret)
        
    } else {
        
        # Univariate potential analysis for all taxa with full data
        potential.results <- list()
        nmodes <- c()
    
        if (is.null(rownames(x))) {
            rownames(x) <- as.character(seq_len(nrow(x)))
        }
        
        for (tax in rownames(x)) {
            if (verbose) {
                message(tax)
            }
            m <- multimodality(as.numeric(x[tax, ]), peak.threshold,
                bw.adjust, bs.iter, min.density, verbose)
            nmodes[[tax]] <- m$modes
            potential.results[[tax]] <- m
        }
        
        multimodality.score <- vapply(potential.results, function(x) {
            1 - x$unimodality.support
        }, 1)
        
        ret <- list(score=multimodality.score,
                    modes=nmodes,
                    results=potential.results)
        
    }
    
    ret
    
}



#' @title Sarle's Bimodality Coefficient
#' @description Sarle's bimodality coefficient.
#' @param x Data vector for which bimodality will be quantified
#' @param bs.iter Bootstrap iterations
#' @param type Score type ('Sarle.finite.sample' or 'Sarle.asymptotic')
#' @return Bimodality score
#' @examples
#' # b <- bimodality_sarle(rnorm(50), type='Sarle.finite.sample')
#' @details The coefficient lies in (0, 1).
#' 
#' The 'Sarle.asymptotic' version is defined as
#' \deqn{b=(g^2 + 1) / k}.
#' This is coefficient of bimodality from Ellison AM Am. J. Bot. 1987, 
#' for microbiome analysis it has been used for instance in
#' Shade et al. 2014.
#'
#' The formula for 'Sarle.finite.sample' (SAS 2012):
#'
#' \deqn{b=\frac{g^2 + 1}{k + (3(n-1)^2)/((n-2)(n-3))}}
#' where n is sample size and 
#' 
#' In both formulas, \eqn{g} is sample skewness and \eqn{k} is the kth
#' standardized moment (also called the sample kurtosis, or
#' excess kurtosis).
#'
#' @references
#' \itemize{
#'   \item{}{Shade et al. mBio 5(4):e01371-14, 2014.}
#'   \item{}{Ellison AM (1987) Am J Botany 74(8):1280-1288.}
#'   \item{}{SAS Institute Inc. (2012). SAS/STAT 12.1 user's guide. Cary, NC.}
#'   \item{}{To cite the microbiome R package, see citation('microbiome')}
#' }
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @seealso Check the dip.test from the \pkg{DIP} package for a
#' classical test of multimodality.
#' @keywords utilities
bimodality_sarle <- function(x, bs.iter=1, type="Sarle.finite.sample") {
    
    g <- skew(x)
    k <- kurtosis(x)
    
    if (type == "Sarle.asymptotic") {
        
        s <- (1 + g^2)/(k + 3)
        
    } else if (type == "Sarle.finite.sample") {
        
        n <- length(x)
        s <- (g^2 + 1)/(k + (3 * (n - 1)^2)/((n - 2) * (n - 3)))
        
    }
    
    if (bs.iter > 1) {
        s <- c()
        for (i in seq_len(bs.iter)) {
            xbs <- sample(x, replace=TRUE)
            s[[i]] <- bimodality_sarle(xbs, type=type)
        }
        s <- mean(s)
    }
    
    s
    
}


# Inspired by moments::kurtosis but rewritten. Internal.
kurtosis <- function (x, na.rm=TRUE)
{
    if (is.matrix(x)) {
        k <- apply(x, 2, kurtosis, na.rm=na.rm)
        return(k)
    } else if (is.vector(x)) {
        if (na.rm) {
            x <- x[!is.na(x)]
    }
        n <- length(x)
        k <- n * sum((x - mean(x))^4)/(sum((x - mean(x))^2)^2)
        return(k)
    } else if (is.data.frame(x)) {
        k <- vapply(x, kurtosis, 1, na.rm = na.rm)
        return(k)
    } else {
        kurtosis(as.vector(x), na.rm = na.rm)
        return(k)
    }
}
