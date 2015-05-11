#' intermediate_stability
#'
#' Quantify stability with respect to a given reference point for all variables.
#'
#' @param dat Input data matrix (variables x samples)
#' @param meta Metadata (samples x factors). This should contain for each sample 
#'           the following self-explanatory fields: subjectID and time. 
#' @param reference.point Calculate stability of the data w.r.t. this point. 
#'                    By default the intermediate range is used (min + (max - min)/2)
#' @param method "lm" (linear model) or "correlation"; the linear model takes time into account 
#' 	         as a covariate 
#'
#' @return A list with following elements: 
#' 	     stability: vector listing the intermediate stability for each variable
#'	     data: processed data sets used for calculations	    
#'
#'
#' @details This method decomposes the data set into differences between
#' consecutive time points. For each variable and time point we calculate
#' for the data values: (i) the distance from reference point; (ii)
#' distance from the data value at the consecutive time point. The
#' "correlation" method calculates correlation between these two
#' variables. Negative correlations indicate that values closer to
#' reference point tend to have larger shifts in the consecutive time
#' point. The "lm" method takes the time lag between the consecutive time
#' points into account as this may affect the comparison and is not taken
#' into account by the straightforward correlation. Here the coefficients
#' of the following linear model are used to assess stability:
#' abs(change) ~ time + abs(start.reference.distance)
#'
#' @export
#'
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples 
#'   # Create simulated example data
#'   dat <- matrix(rnorm(1000), nrow = 10)
#'   rownames(dat) <- paste("Variable", 1:nrow(dat), sep = "")
#'   colnames(dat) <- paste("Sample", 1:ncol(dat), sep = "")
#'   meta <- data.frame(list(
#'     	  sampleID = colnames(dat), 
#'   	  subjectID = rep(paste("subject", 1:50, sep = "-"), each = 2), 
#'	  time = rep(1:2, 50)))
#'   # Intentionally make point very unstable around 0
#'   dat[1, meta$time == 2] <- 1/abs(dat[1, meta$time == 1])
#'   s <- intermediate_stability(dat, meta, reference.point = 0)
#'
intermediate_stability <- function (dat, meta, reference.point = NULL, method = "lm") {

  if (!all(c("subjectID", "time") %in% names(meta))) {
    stop("No subjectID and/or time field provided in metadata!")
  }

  df <- meta
  stabilities <- c()	      
  stabilities.right <- c()	      
  stabilities.left <- c()	      
  data <- list()
  for (i in 1:nrow(dat)) {	      
    df$data <- dat[i,]
    stab <- estimate_stability_single(df = df, reference.point = reference.point, method = method)
    stabilities[[i]] <- stab$stability
    stabilities.right[[i]] <- stab$stability.right
    stabilities.left[[i]] <- stab$stability.left
    data[[i]] <- stab$data
  }

  if (!is.null(rownames(dat))){
    names(stabilities) <- rownames(dat)
    names(stabilities.right) <- rownames(dat)
    names(stabilities.left) <- rownames(dat)
    names(data) <- rownames(dat)
  }

  list(stability = stabilities, stability.right = stabilities.right, stability.left = stabilities.left, data = data)

}

