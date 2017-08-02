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
#' @param phylo.class Whether to return an object of class \code{"phylo"}, in which case a single tree from the 
#'   posterior is returned (not with summary infection times).
#' @return A \code{data.frame} with per item (=node) its parent and support per clade, and optionally summary node times. 
#'   The first n nodes are the samples, the last n-1 nodes are the internal nodes.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' phylotree(MCMCstate)
#' plot(phylotree(MCMCstate, phylo.class = TRUE))
#' @export
phylotree <- function(phybreak.object, samplesize = Inf, support = c("proportion", "count"), phylo.class = FALSE) {
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

    ### initialization
    samplesize <- min(samplesize, chainlength)
    nsamples <- phybreak.object$d$nsamples
    obs <- phybreak.object$p$obs
    res <- c()
    
    ### decision tree to use correct model
    # if(method[1] == 'mcc') {
    res <- .mcctree(.makephyloparsets(phybreak.object$s$nodeparents[, (1:samplesize) + chainlength - samplesize], nsamples), 
                    phybreak.object$s$nodetimes[1:(nsamples - 1), (1:samplesize) + chainlength - samplesize], c(nsamples, samplesize))
    if (phylo.class) 
        return(get.phylo(phybreak.object, samplenr = tail(res, 1) + chainlength - samplesize, simmap = TRUE))
    res <- matrix(head(res, -1), ncol = 5)
    # } if(method[1] == 'cc.construct') { res <- matrix(.CCphylotreeconstruct( .makephyloparsets(phybreak.object$s$nodeparents[,
    # (1:samplesize) + chainlength - samplesize]), phybreak.object$s$nodetimes[1:(obs - 1), (1:samplesize) + chainlength -
    # samplesize], c(obs, samplesize) ), ncol = 4) res[1:obs, 3] <- phybreak.object$v$nodetimes[1:obs] } if(length(res) == 0) {
    # stop('incorrect method provided, choose \'mcc\' or \'cc.construct\'') }
    
    ### make node obs+1 the root node
    if(res[nsamples + 1, 1] != 0) {
      curroot <- which(res[, 1] == 0)
      rootstats <- res[curroot, ]
      res[curroot, ] <- res[nsamples + 1, ]
      res[nsamples + 1, ] <- rootstats
      
      res[res[, 1] == nsamples + 1, 1] <- -1
      res[res[, 1] == curroot, 1] <- nsamples + 1
      res[res[, 1] == -1, 1] <- curroot
    }
    
    # support
    if (support[1] == "count") {
      support.out <- res[, 2]
    } else {
      support.out <- res[, 2]/samplesize
    }
    
    # times
    res[1:nsamples, 3] <- phybreak.object$v$nodetimes[1:nsamples]
    res[1:nsamples, 5] <- phybreak.object$v$nodetimes[1:nsamples]
    
    parents.out <- matrix(res[, 1], ncol = 1, dimnames = list(1:(2 * nsamples - 1), "parent"))
    # if(method[1] == 'mcc') {
    return(data.frame(parents = parents.out, support = support.out, node.times.mean = res[, 3], node.times.sd = res[, 4], node.times.mc.tree = res[, 
        5]))
    # } else { return( data.frame( parents = parents.out, support = res[,2], nodetime.mean = res[,3], nodetime.sd = res[,4] ) ) }
    
    
}


### function to remove the transmission nodes from a matrix with nodeparents-sets called by: .phylotree calls: .makephyloparset
.makephyloparsets <- function(nodeparentsets, nsamples) {
    apply(nodeparentsets, MARGIN = 2, .makephyloparset, nsamples = nsamples)
}

### function to remove the transmission nodes from a nodeparents-set of a phybreak-object called by: .makephyloparsets
.makephyloparset <- function(parentset, nsamples) {
    res <- parentset
    while (max(res) >= 2 * nsamples) {
        res[res >= 2 * nsamples] <- parentset[res][res >= 2 * nsamples]
    }
    return(res[1:(2 * nsamples - 1)])
}

