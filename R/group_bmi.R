#' @title Body-Mass Index (BMI) Classes
#' @description Cut BMI information to standard discrete factors.
#' @param x Numeric vector (BMI)
#' @param breaks Class break points. Either a vector of breakpoints, or one of
#'    the predefined options ("standard", "standard_truncated", "even").
#' @param n Number of groups for the breaks = "even" option.
#' @inheritParams base::cut
#' @return Factor of BMI groups.
#' @references See citation('microbiome') 
#' @author Contact: Leo Lahti \email{microbiome-admin@@googlegroups.com}
#' @export
#' @details Regarding the breaks arguments, the "even" option aims to
#' cut the samples in groups with approximately the same size (by
#' quantiles). The "standard" option corresponds to standard obesity
#' categories defined by the cutoffs <18.5 (underweight); <25 (lean);
#' <30 (obese); <35 (severe obese); <40 (morbid obese); <45 (super
#' obese). The standard_truncated combines the severe, morbid and super
#' obese into a single group.
#' @seealso base::cut
#' @examples
#' bmi.numeric <- range(rnorm(100, mean = 25, sd = 3))
#' bmi.factor <- group_bmi(bmi.numeric)
#' @keywords utilities
group_bmi <- function(x, breaks = "standard", n = 10, labels = NULL,
        include.lowest = TRUE, right = FALSE, dig.lab = 3,
        ordered_result = FALSE) {

    # Predefined groups.
    # Each interval is taken as semi-open to the right by default
    if (length(breaks) == 1 && is.character(breaks)) {    
        if (breaks == "standard") {
            breaks <- c(0, 18.5, 25, 30, 35, 40, 45, Inf)
            labs <- c("underweight",
                "lean", "overweight", "obese", "severe", "morbid", "super")
            return(cut(x, breaks,
                labels = labs, include.lowest = TRUE, right = FALSE))
            # bmi.group <- factor(bmi.group, levels = labs)

        } else if (breaks == "standard_truncated") {
            breaks <- c(0, 18.5, 25, 30, 35, Inf)
            labs <- c("underweight",
                "lean", "overweight", "obese", "severe")
            return(cut(x, breaks,
                labels = labs, include.lowest = TRUE, right = FALSE))
            # bmi.group <- factor(bmi.group, levels = labs)

        } else if (breaks == "even") {
            breaks <- quantile(na.omit(x), seq(0,1,length=n+1))
        } else {
            stop("Input for breaks argument not recognized. 
                See the function help for the available options.")
        }
    }
    
    cut(x, breaks, labels, include.lowest, right, dig.lab, ordered_result)
    
}
