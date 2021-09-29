#' @title gktau
#' @description Measure association between nominal (no order for levels)
#' variables 
#' @details Measure association between nominal (no order for levels) variables 
#' using Goodman and Kruskal tau. Code modified from the original source:
#' r-bloggers.com/measuring-associations-between-non-numeric-variables/
#' An important feature of this procedure is that it allows missing
#' values in either of the variables x or y, treating 'missing' as an
#' additional level.  In practice, this is sometimes very important since
#' missing values in one variable may be strongly associated with either
#' missing values in another variable or specific non-missing levels of
#' that variable. An important characteristic of Goodman and Kruskal's tau 
#' measure is its asymmetry: because the variables x and y enter this expression
#' differently, the value of a(y,x) is not the same as the value of
#' a(x, y), in general.  This stands in marked contrast to either the
#' product-moment correlation coefficient or the Spearman rank
#' correlation coefficient, which are both symmetric, giving the same
#' association between x and y as that between y and x.  The fundamental
#' reason for the asymmetry of the general class of measures defined
#' above is that they quantify the extent to which the variable x is
#' useful in predicting y, which may be very different than the extent to
#' which the variable y is useful in predicting x.
#'
#' @param x first variable
#' @param y second variable
#'
#' @return Dependency measure
#'
#' @examples 
#' data(peerj32)
#' v1 <- factor(peerj32$microbes[,1])
#' v2 <- factor(peerj32$meta$gender)
#' tc <- gktau(v1, v2)
#' 
#' @export
#' @references 
#' Code modified from the original source:
#' \url{
#' http://r-bloggers.com/measuring-associations-between-non-numeric-variables/
#' }
#' To cite the microbiome R package, see citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @keywords utilities
gktau <- function(x, y) {

    if (!is.factor(x)) {
        stop("No factors in x: only factors can be associated with the given method")
    }

    if (!is.factor(y)) {
        stop("No factors in y: only factors can be associated with the given method")
    }
    
    # First, compute the IxJ contingency table between x and y
    Nij <- table(x, y, useNA="ifany")
    
    # Next, convert this table into a joint probability estimate
    PIij <- Nij/sum(Nij)
    
    # Compute the marginal probability estimates
    PIiPlus <- apply(PIij, MARGIN=1, sum)
    PIPlusj <- apply(PIij, MARGIN=2, sum)
    
    # Compute the marginal variation of y
    Vy <- 1 - sum(PIPlusj^2)
    
    # Compute the expected conditional variation of y given x
    InnerSum <- apply(PIij^2, MARGIN=1, sum)
    VyBarx <- 1 - sum(InnerSum/PIiPlus)
    
    # Compute and return Goodman and Kruskal's tau measure
    tau <- (Vy - VyBarx)/Vy
    tau
}
