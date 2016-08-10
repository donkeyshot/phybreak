### consensus phylogenetic tree from phybreak-object ###


### maximum clade credibility tree calls: .mcctree ##C++ .makephyloparsets .makephylo2

#' Maximum clade credibility tree.
#' 
#' Identify the maximum clade credibility tree from a \code{phybreak}-object containing posterior samples.
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param samplesize The number of posterior samples that is used, taken from the tail.
#' @param support Whether to return the support (= posterior probability) for each infector as a \code{"proportion"} 
#'   or as a \code{"count"} of posterior trees in which that transmission link or transmission cluster is present.
#' @param clade.times Whether to include summary times for root nodes of clades, based only on posterior trees 
#'   with the same clade as in the consensus tree. These are used to calculate a mean and standard deviation.
#' @param time.quantiles Used only if \code{clade.times = TRUE}.
#' @param phylo.class Whether to return an object of class \code{"phylo"}, in which case a single tree from the 
#'   posterior is returned (not with summary infection times).
#' @return A \code{data.frame} with per item (=node) its parent and support per clade, and optionally summary node times. 
#'   The first n nodes are the samples, the last n-1 nodes are the internal nodes.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 20)
#' MCMCstate <- phybreak(simulation)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 200)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 100, thin = 10)
#' 
#' phylotree(MCMCstate)
#' @export
phylotree <- function(phybreak.object, samplesize = Inf, support = c("proportion", "count"), clade.times = TRUE, time.quantiles = c(0.025, 
    0.5, 0.975), phylo.class = FALSE) {
    chainlength <- length(phybreak.object$s$mu)
    
    ### tests
    if (chainlength == 0) 
        stop("no sampled trees available")
    if (samplesize > chainlength & samplesize < Inf) {
        warning("desired 'samplesize' larger than number of available samples")
    }
    if (support[1] != "proportion" && support[1] != "count") {
        warning("support is given as proportion")
    }
    if (class(time.quantiles) != "numeric") {
        stop("time.quantiles should be numeric")
    }
    if (max(time.quantiles) > 1 || min(time.quantiles) < 0) {
        stop("time.quantiles should be between 0 and 1")
    }
    
    samplesize <- min(samplesize, chainlength)
    obs <- phybreak.object$p$obs
    
    res <- c()
    
    # if(method[1] == 'mcc') {
    res <- .mcctree(.makephyloparsets(phybreak.object$s$nodeparents[, (1:samplesize) + chainlength - samplesize]), phybreak.object$s$nodetimes[1:(obs - 
        1), (1:samplesize) + chainlength - samplesize], c(obs, samplesize))
    if (phylo.class) 
        return(get.phylo(phybreak.object, tail(res, 1) + chainlength - samplesize, TRUE))
    res <- matrix(head(res, -1), ncol = 5)
    res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs]
    res[1:obs, 5] <- phybreak.object$v$nodetimes[1:obs]
    # } if(method[1] == 'cc.construct') { res <- matrix(.CCphylotreeconstruct( .makephyloparsets(phybreak.object$s$nodeparents[,
    # (1:samplesize) + chainlength - samplesize]), phybreak.object$s$nodetimes[1:(obs - 1), (1:samplesize) + chainlength -
    # samplesize], c(obs, samplesize) ), ncol = 4) res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs] } if(length(res) == 0) {
    # stop('incorrect method provided, choose \'mcc\' or \'cc.construct\'') }
    
    parents.out <- matrix(res[, 1], ncol = 1, dimnames = list(1:(2 * obs - 1), "parent"))
    # if(method[1] == 'mcc') {
    return(data.frame(parents = parents.out, support = res[, 2], node.times.mean = res[, 3], node.times.sd = res[, 4], node.times.mc.tree = res[, 
        5]))
    # } else { return( data.frame( parents = parents.out, support = res[,2], nodetime.mean = res[,3], nodetime.sd = res[,4] ) ) }
    
    
}


### function to remove the transmission nodes from a matrix with nodeparents-sets called by: .phylotree calls: .makephyloparset
.makephyloparsets <- function(nodeparentsets) {
    apply(nodeparentsets, MARGIN = 2, .makephyloparset)
}

### function to remove the transmission nodes from a nodeparents-set of a phybreak-object called by: .makephyloparsets
.makephyloparset <- function(parentset) {
    obs <- (1 + length(parentset))/3
    res <- parentset
    while (max(res) >= 2 * obs) {
        res[res >= 2 * obs] <- parentset[res][res >= 2 * obs]
    }
    return(res[1:(2 * obs - 1)])
}

