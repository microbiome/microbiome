#' Description: asses significance of hierarchical clustering 
#' clusters based on multiscale bootstrap resamping
#' 
#' @param dat data matrix samples to be clustered on the rows
#' @param my.groups groups from hierarchical clustering
#' @param R bootstrap iterations
#' @param sample.sizes sample sizes for multiscale bootstrap as fractions of 
#'                  the total data (0...1..etc.). 
#'               Use vector of different sizes to perform 
#'               multiscale bootstrap resampling.
#' @param min.size minimum cluster size
#' @param corr.th correlation threshold for determining hclust groups
#' @param replace sampling with replacement TRUE/FALSE
#' @param metric correlation metrix
#' @param verbose verbose
#' @param pvalue.threshold pvalue threshold: filter out non-significant groups
#' @param remove.nested.clusters Remove nested clusters, or keep all? 
#'                    Default: TRUE
#'
#' @return List of significant clusters (corresponding to the given p-value)
#'
#' @export
#'
#' @references See citation('microbiome') 
#' @author Leo Lahti \email{leo.lahti@@iki.fi}
#' @examples data(peerj32); 
#'          sign.clust <- hclust.significance(
#'                        as.matrix(peerj32$microbes), 
#'                   R = 20, 
#'                   min.size = 2, 
#'                   corr.th = 0.5) 
#' @keywords methods

hclust.significance <- function(dat, my.groups = NULL, R, sample.sizes = 1, 
                                min.size, 
    corr.th, replace = TRUE, metric = "pearson", verbose = TRUE, 
    pvalue.threshold = 0.05, 
    remove.nested.clusters = TRUE) {
    
    if (is.null(my.groups)) {
        if (verbose) {
            message("Forming the clusters")
        }
        my.groups <- get.hclust.groups(dat, corr.th = corr.th, 
                          recursive = TRUE, 
            min.size = min.size, metric = metric)
    }
    
    # Group the phylotypes based on random subsets of the data
    checks <- matrix(NA, ncol = length(my.groups), 
                nrow = R * length(sample.sizes))
    
    if (verbose) {
        message("Bootstrap resampling")
    }
    # Use multiscale bootstrap resampling to get more robust p-value estimates
    cnt <- 0
    for (sample.size in round(sample.sizes * ncol(dat))) {
        
        for (r in 1:R) {
            
            if (verbose) {
                message(paste(sample.size, " ", r))
            }
            
            # Bootstrap sampling
            inds <- sample(ncol(dat), replace = replace, size = sample.size)
            dat2 <- dat[, inds]
            # Ensure std is never zero for items as a result of resampling
            while (any(apply(dat2, 1, sd) == 0)) {
                inds <- sample(ncol(dat), replace = replace, size = sample.size)
                dat2 <- dat[, inds]
            }
            
            # Get groups that exceed threshold. The most general
            # groups, no nested ones.

            my.groups2 <- get.hclust.groups(dat2, corr.th, recursive = FALSE, 
                             min.size = min.size)
            
            # For each group in the original data, check if it is
            # included in one of the bootstrap groups

            check <- c()
            for (i in 1:length(my.groups)) {
                g1 <- my.groups[[i]]
                verified <- FALSE
                for (g2 in my.groups2) {
                  if (all(g1 %in% g2)) {
                    verified <- TRUE
                    break
                  }
                }
                check[[i]] <- verified
            }
            
            cnt <- cnt + 1
            checks[cnt, ] <- check
            
        }
    }
    
    pvals <- 1 - colMeans(checks)
    
    if (verbose) {
        message("Correct p-values for multiple testing")
    }
    pvals <- p.adjust(pvals, "BH")
    
    # Pick element names instead of indices
    my.groups <- lapply(my.groups, function(inds) {
        rownames(dat)[inds]
    })
    
    # Pick the significant clusters
    significant.clusters <- my.groups[which(pvals < pvalue.threshold)]
    
    if (remove.nested.clusters) {
        if (verbose) {
            message("Removing nested clusters")
        }
        significant.clusters <- remove.nested.groups(significant.clusters)
    }
    
    significant.clusters
    
}


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
