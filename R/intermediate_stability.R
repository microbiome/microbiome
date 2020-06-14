#' @title Intermediate Stability
#' @description Quantify intermediate stability with respect to a given
#' reference point. 
#' @param x \pkg{phyloseq} object.
#' Includes abundances (variables x samples) and
#' sample_data data.frame (samples x features) with 'subject'
#' and 'time' field for each sample.
#' @param reference.point Calculate stability of the data w.r.t. this point.
#' By default the intermediate range is used (min + (max - min)/2).
#' If a vector of points is provided, then the scores will be calculated
#' for every point and a data.frame is returned.
#' @param method 'lm' (linear model) or 'correlation';
#' the linear model takes time into account as a covariate 
#' @param output Specify the return mode. Either the 'full' set of stability
#' analysis outputs, or the 'scores' of intermediate stability.
#' @return A list with following elements: 
#' stability: estimated stability
#' data: processed data set used in calculations        
#' @details Decomposes each column in x into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' 'correlation' method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The 'lm' method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance). Samples with missing
#' data, and subjects with less than two time point are excluded. The absolute
#' count data x is logarithmized before the analysis with the log10(1 + x)
#' trick to circumvent logarithmization of zeroes.
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @export
#' @examples
#' data(atlas1006)
#' x <- subset_samples(atlas1006, DNA_extraction_method == 'r')
#' x <- prune_taxa(c('Akkermansia', 'Dialister'), x)
#' res <- intermediate_stability(x, reference.point=NULL)
#' @keywords utilities
intermediate_stability <- function(x, reference.point=NULL,
    method="correlation", output="scores") {
    
    x0 <- x
    
    if (length(reference.point) > 1 && output == "scores") {
        
        scores <- c()
        for (i in seq_len(length(reference.point))) {
            scores[[i]] <- intermediate_stability(x,
        reference.point=reference.point[[i]], 
                method=method, output=output)
        }
        
        if (ntaxa(x) > 1) {
            scores <- as.data.frame(scores, nrow=ntaxa(x))
        } else {
            scores <- as.data.frame(t(as.matrix(scores,
            ncol=length(reference.point))))
        }
        
        colnames(scores) <- as.character(reference.point)
        rownames(scores) <- taxa(x)
        
        return(scores)
    }
    
    # Logarithmize the data with CLR 
    x <- abundances(transform(x0, "clr"))
    meta <- sample_data(x0)
    
    # Estimate stabilities for each OTU
    stability <- list()
    df <- meta
    
    # Ensure time is numeric
    df$time <- as.numeric(as.character(df$time))
    
    # Remove subjects with only one measurement
    df <- df[df$subject %in% names(which(table(df$subject) > 1)), ]
    
    # Split data by subject
    spl <- split(df, as.character(df$subject))
    spl.list <- list()
    for (i in seq_len(length(spl))) {
        # Ensure the measurements are ordered in time
        spl.list[[i]] <- list(spl=spl[[i]][order(spl[[i]]$time), ],
        time.difs=diff(spl[[i]]$time))
    }
    
    for (tax in rownames(x)) {
        
        df$data <- x[tax, rownames(df)]
        # Remove NAs and Infinities keep <- which(!is.na(df$data))
        # df <- df[keep,]
        
        stability[[tax]] <- estimate_stability(df,
        reference.point=reference.point, 
            method=method, spl.list)
        
    }
    
    if (output == "full") {
        return(stability)
    } else if (output == "scores") {
        intermediate.stability <- vapply(stability, function(x) {
            x$stability
        }, 1)
        return(intermediate.stability)
    }
    
}




#' @title Estimate Stability
#' @description Quantify intermediate stability with respect to a given
#'    reference point. 
#' @param df Combined input data vector (samples x variables) and metadata
#'        data.frame (samples x features)
#'        with the 'data', 'subject' and 'time' field for each sample 
#' @param reference.point Optional. Calculate stability of the data w.r.t.
#'                        this point. By default the intermediate range is
#'                        used (min + (max - min)/2)
#' @param method 'lm' (linear model) or 'correlation'; the linear model takes
#'            time into account as a covariate
#' @param spl split object to speed up
#' @return A list with following elements: 
#'        stability: estimated stability
#'        data: processed data set used in calculations        
#' @details Decomposes each column in x into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' 'correlation' method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The 'lm' method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance). 
#' Samples with missing data, and subjects with less than two time point are
#' excluded.       
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples
#' # df <- data.frame(list(
#' #           subject=rep(paste('subject', 1:50, sep='-'), each=2), 
#' #           time=rep(1:2, 50),
#' #           data=rnorm(100)))
#' #s <- estimate_stability_single(df, reference.point=NULL, method='lm')
#' @keywords internal
estimate_stability <- function(df, reference.point=NULL, method="lm",
    spl.list) {
    
    # Detect intermediate value in the overall data if reference point
    # not given
    if (is.null(reference.point)) {
        reference.point <- mean(range(df$data))
    }
    
    if (nrow(df) < 2) {
        warning("No subjects with time series in estimate_stability. 
        Returninng NULL")
        return(NULL)
    }
    
    dfis <- NULL
    
    for (i in seq_len(length(spl.list))) {
        
        spli <- spl.list[[i]]$spl
        spli$data <- unlist(df[rownames(spli), "data"], use.names=FALSE)
        time.difs <- spl.list[[i]]$time.difs
        
        # Calculate differences between consecutive time points;
        # and the initial values;
        # and their distance from reference
        data.difs <- diff(spli$data)
        
        start.points <- spli$data[-nrow(spli)]
        start.reference.distance <- start.points - reference.point
        
        # Organize into data frame
        dfi <- data.frame(change=data.difs,
        time=time.difs, start=start.points, 
            start.reference.distance=start.reference.distance)
        
        # Add to the collection
        dfis <- rbind(dfis, dfi)
        
    }
    
    dfis.left <- subset(dfis, start.reference.distance < 0)
    dfis.right <- subset(dfis, start.reference.distance > 0)
    
    # Simplified stability calculation (do not consider time effect)
    stability <- NULL
    stability.left <- stability.right <- NA
    
    if (method == "correlation") {
        
        # For each subject, check distance from the stability point at
        # the baseline time point
    
        baseline.distance <- abs(dfis$start.reference.distance)
        
        # For each subject, calculate deviation between the first and
        # second time point
    
        followup.distance <- abs(dfis$change)
        stability <- cor(baseline.distance, followup.distance)
        
        if (nrow(dfis.left) > 10) {
            # Negative values for low stability
            baseline.distance <- abs(dfis.left$start.reference.distance)
            followup.distance <- dfis.left$change
            stability.left <- cor(baseline.distance, followup.distance)
        }
        
        if (nrow(dfis.right) > 10) {
            baseline.distance <- abs(dfis.right$start.reference.distance)
            followup.distance <- dfis.right$change
            stability.right <- -cor(baseline.distance, followup.distance)
        }
        
        
    } else if (method == "lm") {
    
        # Advanced calculation, take time into account with linear
        # model (also possible to check p-values later if needed)
    
        stability <- coef(summary(lm(abs(change) ~ time +
        abs(start.reference.distance), 
            data=dfis)))["abs(start.reference.distance)", "Estimate"]
        
        if (nrow(dfis.left) > 10) {
            # Negative values for low stability
            stability.left <- coef(summary(lm(change ~ time +
            abs(start.reference.distance), 
                data=dfis.left)))["abs(start.reference.distance)", "Estimate"]
        }
        if (nrow(dfis.right) > 10) {
            # Negative values for low stability
            stability.right <- -coef(summary(lm(change ~ time +
            abs(start.reference.distance), 
            data=dfis.right)))["abs(start.reference.distance)", "Estimate"]
            
        }
    }
    
    list(stability=stability,
        stability.right=stability.right, stability.left=stability.left, 
        data=dfis)
    
}
