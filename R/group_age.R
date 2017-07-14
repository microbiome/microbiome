#' @title Age Classes
#' @description Cut age information to discrete factors.
#' @param x Numeric vector (age in years)
#' @param breaks Class break points. Either a vector of breakpoints,
#' or one of the predefined options ("years", "decades", "even").
#' @param n Number of groups for the breaks = "even" option.
#' @inheritParams base::cut
#' @return Factor of age groups.
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @details Regarding the breaks arguments, the "even" option aims to
#' cut the samples in groups with approximately the same size (by
#' quantiles). The "years" and "decades" options are self-explanatory.
#' @seealso base::cut
#' @examples
#' data(atlas1006)
#' age.numeric <- meta(atlas1006)$age
#' age.factor <- group_age(age.numeric)
#' @keywords utilities
group_age <- function(x, breaks = "decades", n = 10, labels = NULL,
        include.lowest = TRUE, right = FALSE, dig.lab = 3,
        ordered_result = FALSE) {

    # Predefined groups.
    # Each interval is taken as semi-open to the right by default
    if (length(breaks) == 1 && is.character(breaks)) {
        if (breaks == "decades") {
            breaks <- seq(floor(min(na.omit(x))/10)*10,
                    ceiling(max(na.omit(x))/10)*10, 10)
        } else if (breaks == "years") {
            breaks <- seq(floor(min(na.omit(x))),
                    ceiling(max(na.omit(x))), 1)
        } else if (breaks == "even") {
            breaks <- quantile(na.omit(x), seq(0,1,length=n+1))
        }
    }
    
    cut(x, breaks, labels, include.lowest, right, dig.lab, ordered_result)

}
