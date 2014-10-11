
#' Description: remove nested groups from hierarchical clustering results; 
#'              keep only the largest group when nested groups are detected
#' 
#'
#' @param groups list of groups
#'
#' @return list of groups with nested groups filtered out
#'
#' @export
#'
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples data(peerj32); 
#'          cl <- get.hclust.groups(peerj32$microbes, 
#'                                 corr.th = 0.8)
#'          cl <- remove.nested.groups(cl)
#' @keywords methods

remove.nested.groups <- function(groups) {
    
    # Remove nested groups
    nested <- rep(FALSE, length(groups))
    for (i in 1:length(groups)) {
        for (j in setdiff(1:length(groups), i)) {
            if (all(groups[[i]] %in% groups[[j]])) {
                nested[[i]] <- TRUE
            }
        }
    }
    groups[!nested]
    
}


#' Description: Get hierarchical clustering groups, 
#' including the nested ones if recursive = TRUE
#' 
#'
#' @param dat data matrix
#' @param corr.th correlation threshold for hclust
#' @param recursive get also nested groups TRUE/FALSE
#' @param min.size minimum cluster size
#' @param metric similarity measure (spearman / pearson)
#'
#' @return groups
#'
#' @import fastcluster
#'
#' @export
#'
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples data(peerj32); 
#'          cl <- get.hclust.groups(peerj32$microbes, 
#'                                 corr.th = 0.8)
#' @keywords methods

get.hclust.groups <- function(dat, corr.th, recursive = FALSE, min.size = 2, 
                              metric = "pearson") {
    
    # Group the phylotypes
    cordist <- as.dist(1 - cor(t(dat), method = metric))
    hc <- hclust(cordist, method = "complete")
    
    # Identify all groups that are higher than threshold First,
    # identify number of groups at the threshold

    groups <- cutree(hc, h = 1 - corr.th)
    K <- length(unique(groups))
    
    # Then, get all nested clusterings (ie those with more groups)
    if (recursive) {
        
        my.groups <- list()
        cnt <- 0
        for (k in K:(nrow(dat) - 1)) {
            groups <- cutree(hc, k = k)
            
            # Pick groups with > 1 members
            uniq.groups <- unique(groups)[sapply(unique(groups), function(g) {
                sum(groups == g)
            }) >= min.size]
            for (i in 1:length(uniq.groups)) {
                cnt <- cnt + 1
                g <- uniq.groups[[i]]
                my.groups[[cnt]] <- which(groups %in% g)
            }
        }
        
        # Remove duplicated groups
        my.groups <- my.groups[!duplicated(my.groups)]
        
    } else {
        
        my.groups <- list()
        uniq.groups <- unique(groups)[sapply(unique(groups), function(g) {
            sum(groups == g)
        }) >= min.size]
        for (i in 1:length(uniq.groups)) {
            g <- uniq.groups[[i]]
            my.groups[[i]] <- which(groups %in% g)
        }
        
    }
    
    my.groups
} 
