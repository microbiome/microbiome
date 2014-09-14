#' Description: ROC analysis returning true and false positive rates
#' along an ordered list
#'
#' Arguments:
#' @param ordered.results Items ordered from best to worst according
#'        to the test score.
#' @param true.positives known true positives
#'
#' Returns:
#'   @return List: true positive rate (tpr) and false positive rate (fpr)
#'
#' @export
#'
#' @examples data(peerj32); 
#'           x <- unlist(peerj32$microbes[1,]); 
#'	     res <- roc(names(x), sample(names(x), 10))
#'
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

roc <- function (ordered.results, true.positives) {
	
	# Check that all known positives are included in the original 
	# analysis i.e. ordered results list
	positives<-true.positives[true.positives %in% ordered.results]	
	
	# Number of retrieved known cytobands
	N <- length(ordered.results) #total number of samples
	Np <- length(positives) #number of positives
	Nn <- N - Np #number of negatives

	TP <- cumsum(ordered.results %in% positives)
	FP <- cumsum(!(ordered.results %in% positives))
	tpr <- TP/Np #TP/(TP + FN) = TP.simCCA / true.positives
	fpr <- FP/Nn #FP/(FP + TN) = FP.simCCA / N.simCCA

	list(tpr=tpr,fpr=fpr)
}


#' Description: Plot ROC curve
#'
#' Arguments:
#' @param ordered.results Items ordered from best to worst according
#'        to the test score.
#' @param true.positives known true positives
#' @param line Draw 45 angle line
#' @param title Title text
#' Returns:
#'   @return Used for its side effects (plot)
#'
#' @export
#'
#' @examples data(peerj32); 
#' 	     x <- unlist(peerj32$microbes[1,]); 
#' 	     res <- roc.plot(names(x), sample(names(x), 10)) 
#' 
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
roc.plot <- function(ordered.results, true.positives, line=F, title="") {
	res <- roc(ordered.results, true.positives)
	plot(res$fpr,res$tpr,lty=1,type='l',xlab="False positive rate",ylab="True positive rate",xlim=c(0,1),ylim=c(0,1),main=paste("ROC curve", title))
	if (line) {
	  # Draw 45 angle line
	  lines(c(0,1),c(0,1))
	}
}


#' Description: ROC AUC calculation
#'
#' Arguments:
#' @param ordered.results Items ordered from best to worst according
#'        to the test score.
#' @param true.positives known true positives
#'
#' Returns:
#'   @return ROC AUC value
#'
#' @export
#'
#' @examples data(peerj32); 
#' 	     x <- unlist(peerj32$microbes[1,]); 
#'	     res <- roc.auc(names(x), sample(names(x), 10)) 
#' 
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

roc.auc <- function (ordered.results, true.positives) {

  # Compute area under curve
  rates <- roc(ordered.results, true.positives)
  # integration: Compute step intervals and compute weighted sum of 
  # true positive rates in each interval.
  # note that interval length can be 0 if fpr does not change
  auc <- as.numeric(t(rates$fpr[-1]-rates$fpr[-length(rates$fpr)])%*%rates$tpr[-length(rates$fpr)])

  auc
}



#' Description: ROC AUC calculation for a matrix of variables 
#'
#' Arguments:
#' @param dat Data matrix (variables x samples)
#' @param true.positives known true positive samples
#'
#' @details The samples are ordered for each row (variable) from the highest to the lowest score, and ROC/AUC value is calculated based on this ordering.
#'
#' Returns:
#'   @return Vector of ROC AUC values for each variable
#'
#' @export
#'
#' @examples \dontrun{rocs <- roc.auc.matrix(dat, true.positives)}
#' 
#' @references See citation("microbiome") 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities

roc.auc.matrix <- function (dat, true.positives) {

  aucs <- c()
  for (k in 1:nrow(dat)) {

    # Order the scores for this variables from the highest to the lowest
    ordered.results <- names(rev(sort(dat[k, ])))

    # Calculate AUC and store
    auc <- roc.auc(ordered.results, true.positives)
    aucs[[k]] <- auc

  }

  if (!is.null(rownames(dat))) {
    names(aucs) <- rownames(dat)
  }

  aucs

}
