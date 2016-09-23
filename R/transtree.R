### consensus transmission tree from phybreak-object ###


### several methods available, default is 'count' = counting infectors without resolving cycles or multiple/no index cases.
### Other options are 'edmunds' = counting infectors with resolving cycles and multiple/no index cases; 'mpc' = maximum parent
### credibility, selecting the sampled posterior tree with best support (summed count of infectors per host); 'mtcc' = maximum
### transmission cluster credibility, selecting the sampled posterior tree with best support (summed count of clusters, where
### cluster = host + all progeny) calls: .transtreecount .transtreeedmonds .mpcinfector .mtcctree ##C++ get.phylo


#' Create a consensus transmission tree.
#' 
#' Various methods to create summary transmission trees from a \code{phybreak}-object containing posterior samples.
#' 
#' Four methods are supported for transmission tree reconstruction (who infected whom). 
#' \itemize{
#'   \item \code{"count"} gives the most frequent infector from the posterior samples. Note that this may 
#'     result in an improper tree (multiple roots and/or cycles). Support is measured by the frequency of the 
#'     infector in the posterior distribution.
#'   \item \code{"edmonds"} starts from the most frequent infector (method \code{"count"}), multiple roots and 
#'     cycles are removed by selecting one by one the next most frequent option that minimizes the loss in support 
#'     (\href{https://en.wikipedia.org/wiki/Edmonds\%27_algorithm}{Edmonds' algorithm}). Support is measured by the frequency of the infector in the posterior distribution.
#'   \item{"mpc"} gives the maximum parent credibility tree as described in \href{http://journals.plos.org/ploscompbiol/article?id=10.1371/journal.pcbi.1004613}{Hall et al (2015)}. This is the tree 
#'     in the set of posterior samples that has maximum support = product of frequencies among all posterior samples. 
#'     Support is measured by the frequency of the infector in the posterior distribution.
#'   \item{"mtcc"} gives the maximum transmission cluster credibility tree. This is equivalent to the maximum clade
#'     credibility (mcc) phylogenetic tree, with clusters defined as host + all progeny. Support is measured by the 
#'     frequency of the cluster in the posterior distribution.
#' }
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param method The method used to create the tree (see details).
#' @param samplesize The number of posterior samples that is used, taken from the tail.
#' @param infector.name Whether to return the infector names, or the position in the vector of hosts.
#' @param support Whether to return the support (= posterior probability) for each infector as a \code{"proportion"} 
#'   or as a \code{"count"} of posterior trees in which that transmission link or transmission cluster is present.
#' @param infection.times Whether to base the summary infection times on \code{"all"} samples, or only on samples 
#'   in which the \code{"infector"} was the same as in the consensus tree. A third option is \code{"infector.sd"}, 
#'   in which case only posterior trees with the same infector or transmission cluster as in the consensus tree 
#'   were used to calculate a mean and standard deviation. In that case and if \code{method = "mpc"} or \code{method = "mtcc"},
#'   also the infection times in the first sampled tree with consensus tree topology are given.
#' @param time.quantiles Used only if \code{infection.times = "all"} or \code{"infector"}.
#' @param phylo.class Whether to return an object of class \code{"phylo"}, in which case a single tree 
#'   (\code{"mpc"} or \code{"mtcc"}) from the posterior is returned (not with summary infection times).
#' @return If \code{phylo.class = FALSE}, a \code{data.frame} with per item (=host) its infector and support per 
#'   infector (or cluster), and summary infection times. If \code{phylo.class = TRUE}, a class \code{"phylo"} object, a 
#'   single tree (\code{"mpc"} or \code{"mtcc"}) from the posterior is returned (not with summary infection times).
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' transtree(MCMCstate, method = "edmonds")
#' transtree(MCMCstate, method = "mpc", infection.times = "infector.sd")
#' plot(MCMCstate, plot.which = "mpc")
#' @export
transtree <- function(phybreak.object, method = c("count", "edmonds", "mpc", "mtcc"), samplesize = Inf, infector.name = TRUE, 
    support = c("proportion", "count"), infection.times = c("all", "infector", "infector.sd"), time.quantiles = c(0.025, 0.5, 0.975), phylo.class = FALSE) {
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
    if (infection.times[1] != "all" && infection.times[1] != "infector" && infection.times[1] != "infector.sd") {
        warning("infection time summaries based on all samples")
    }
    if (infection.times[1] == "all" || infection.times[1] == "infector") {
        if (class(time.quantiles) != "numeric") {
            stop("time.quantiles should be numeric")
        }
        if (max(time.quantiles) > 1 || min(time.quantiles) < 0) {
            stop("time.quantiles should be between 0 and 1")
        }
    }
    
    ### initialization
    samplesize <- min(samplesize, chainlength)
    obs <- phybreak.object$p$obs
    res <- c()
    
    ### decision tree to use correct model
    if (method[1] == "count") {
        res <- .transtreecount(phybreak.object, samplesize, infection.times[1] == "infector.sd")
    }
    if (method[1] == "edmonds") {
        res <- .transtreeedmonds(phybreak.object, samplesize, infection.times[1] == "infector.sd")
    }
    if (method[1] == "mpc") {
        res <- .mpcinfector(phybreak.object, samplesize, phylo.class, infection.times[1] == "infector.sd")
        if (phylo.class) 
            return(get.phylo(phybreak.object, res, TRUE))
    }
    if (method[1] == "mtcc") {
        res <- .mtcctree(phybreak.object$s$nodehosts[obs:(2 * obs - 1), (1:samplesize) + chainlength - samplesize], phybreak.object$s$nodetimes[obs:(2 * 
            obs - 1), (1:samplesize) + chainlength - samplesize], c(obs, samplesize))
        if (phylo.class) 
            return(get.phylo(phybreak.object, tail(res, 1) + chainlength - samplesize, TRUE))
        res <- matrix(head(res, -1), ncol = 5)
    }
    # if(method[1] == 'cc.construct') { res <- matrix(.CCtranstreeconstruct( phybreak.object$s$nodehosts[obs:(2*obs-1),
    # (1:samplesize) + chainlength - samplesize], phybreak.object$s$nodetimes[obs:(2*obs-1), (1:samplesize) + chainlength -
    # samplesize], c(obs, samplesize) ), ncol = 4) } test: no correct method provided
    if (length(res) == 0) {
        stop("incorrect method provided, choose \"count\", \"edmonds\",
         \"mpc\", or \"mtcc\"")
    }
    
    ### build output infectors
    if (infector.name) {
        infectors.out <- matrix(c("index", phybreak.object$d$names)[1 + res[, 1]], ncol = 1, dimnames = list(phybreak.object$d$names, 
            "infector"))
    } else {
        infectors.out <- matrix(res[, 1], ncol = 1, dimnames = list(phybreak.object$d$names, "infector"))
    }
    # support
    if (support[1] == "count") {
        support.out <- res[, 2]
    } else {
        support.out <- res[, 2]/samplesize
    }
    # times
    if (infection.times[1] == "infector") {
        posttimes <- phybreak.object$s$nodetimes[obs:(2 * obs - 1), (1:samplesize) + chainlength - samplesize]
        posttimes[res[, 1] != phybreak.object$s$nodehosts[obs:(2 * obs - 1), (1:samplesize) + chainlength - samplesize]] <- NA
        time.out <- t(matrix(apply(posttimes, 1, quantile, probs = time.quantiles, na.rm = TRUE), ncol = obs))
        colnames(time.out) <- paste0("inf.times.Q", 100 * time.quantiles)
    } else if (infection.times[1] == "infector.sd") {
        if (method[1] == "mpc" || method[1] == "mtcc") {
            time.out <- res[, 3:5]
            colnames(time.out) <- paste0("inf.times.", c("mean", "sd", "mc.tree"))
        } else {
            time.out <- res[, 3:4]
            colnames(time.out) <- paste0("inf.times.", c("mean", "sd"))
        }
    } else {
        posttimes <- phybreak.object$s$nodetimes[obs:(2 * obs - 1), (1:samplesize) + chainlength - samplesize]
        time.out <- t(matrix(apply(posttimes, 1, quantile, probs = time.quantiles), ncol = obs))
        colnames(time.out) <- paste0("inf.times.Q", 100 * time.quantiles)
    }
    
    ### return the result
    return(data.frame(infectors = infectors.out, support = support.out, time.out))
}


### function to obtain the most likely infector for each host, possibly containing cycles and multiple roots, plus (if
### requested) means and standard deviations of infection times called by: transtree calls: .postinfector
.transtreecount <- function(phybreak.object, samplesize, includetimes) {
    ### initialize some constants
    chainlength <- length(phybreak.object$s$mu)
    obsize <- phybreak.object$p$obs
    samplerange <- (chainlength - samplesize + 1):chainlength
    
    ### get result
    res <- t(matrix(with(phybreak.object, apply(s$nodehosts[obsize:(2 * obsize - 1), samplerange], 1, .postinfector, support = TRUE)), 
        nrow = 2))
    
    ### get time summaries
    if (includetimes) {
        timesums <- with(phybreak.object, rowSums(s$nodetimes[obsize:(2 * obsize - 1), samplerange] * (s$nodehosts[obsize:(2 * 
            obsize - 1), samplerange] == res[, 1])))
        timesumsqs <- with(phybreak.object, rowSums((s$nodetimes[obsize:(2 * obsize - 1), samplerange]^2) * (s$nodehosts[obsize:(2 * 
            obsize - 1), samplerange] == res[, 1])))
        res <- cbind(res, timesums/res[, 2])
        res <- cbind(res, sqrt((timesumsqs - res[, 3]^2 * res[, 2])/(res[, 2] - 1)))
    }
    
    ### return result
    return(res)
}


### function to obtain a transmission tree based on most likely infectors, using Edmonds's algorithm to remove cycles, applied
### with multiple root candidates, plus (if requested) means and standard deviations of infection times called by transtree
### calls: .edmondsiterative
.transtreeedmonds <- function(phybreak.object, samplesize, includetimes) {
    ### initialize some constants
    chainlength <- length(phybreak.object$s$mu)
    obsize <- phybreak.object$p$obs
    samplerange <- (chainlength - samplesize + 1):chainlength
    
    ### obtaining the result in steps
    
    # matrix with support for each infector (row) per host (column), with 0 as maximum, and a column for the index
    supportmatrix <- cbind(c(0, rep(-1, obsize)), apply(1 + phybreak.object$s$nodehosts[obsize:(2 * obsize - 1), samplerange], 
        1, tabulate, nbins = obsize + 1))
    supportmatrix <- supportmatrix - rep(apply(supportmatrix, 2, max), each = obsize + 1)
    
    # vector with hosts, ordered by support to be index, and vector with these supports
    candidateindex <- order(apply(phybreak.object$s$nodehosts[obsize:(2 * obsize - 1), samplerange] == 0, 1, sum), decreasing = TRUE)
    indexsupports <- apply(phybreak.object$s$nodehosts[obsize:(2 * obsize - 1), samplerange] == 0, 1, sum)[candidateindex]
    
    ## index host candidate by candidate, make a tree with edmonds's algorithm after each tree, test for each remaining candidate
    ## if their support-to-be-index is smaller than the support for their infector in the last tree.  If so, stop making new
    ## trees, as no better tree will come out. Then select the best tree.  first initialize some variables
    bestYN <- FALSE
    nextcandidate <- 1
    alltrees <- c()
    allsupports <- c()
    # then make trees as long as bestYN == FALSE
    while (!bestYN) {
        # make copy of supportmatrix, maximally supporting the next candidate index
        suppmat <- supportmatrix
        suppmat[1, ] <- -samplesize
        suppmat[1, candidateindex[nextcandidate] + 1] <- 0
        
        # make the tree
        thistree <- .edmondsiterative(suppmat, samplesize, obsize)[-1] - 1
        alltrees <- c(alltrees, thistree)
        thissupport <- rowSums(phybreak.object$s$nodehosts[obsize:(2 * obsize - 1), samplerange] == thistree)
        allsupports <- c(allsupports, thissupport)
        
        # test if any unused candidate index has higher index support than infector support in last tree
        bestYN <- TRUE
        for (i in tail(candidateindex, -nextcandidate)) {
            if (thissupport[i] < indexsupports[candidateindex == i]) {
                bestYN <- FALSE
            }
        }
        nextcandidate <- nextcandidate + 1
    }
    
    # find the tree with maximum support
    dim(alltrees) <- c(obsize, length(alltrees)/obsize)
    dim(allsupports) <- dim(alltrees)
    besttree <- alltrees[, which(colSums(allsupports) == max(colSums(allsupports)))]
    treesupport <- allsupports[, which(colSums(allsupports) == max(colSums(allsupports)))]
    res <- matrix(c(besttree, treesupport), ncol = 2)
    
    ### get time summaries
    if (includetimes) {
        timesums <- with(phybreak.object, rowSums(s$nodetimes[obsize:(2 * obsize - 1), samplerange] * (s$nodehosts[obsize:(2 * 
            obsize - 1), samplerange] == res[, 1])))
        timesumsqs <- with(phybreak.object, rowSums((s$nodetimes[obsize:(2 * obsize - 1), samplerange]^2) * (s$nodehosts[obsize:(2 * 
            obsize - 1), samplerange] == res[, 1])))
        res <- cbind(res, timesums/res[, 2])
        res <- cbind(res, sqrt((timesumsqs - res[, 3]^2 * res[, 2])/(res[, 2] - 1)))
    }
    
    ### return result
    return(res)
}


### function to obtain a transmission tree based on most likely infectors, selecting the tree among posterior trees with
### highest support, plus (if requested) means and standard deviations of infection times called by transtree
.mpcinfector <- function(phybreak.object, samplesize, phylo.class = FALSE, includetimes = FALSE) {
    ### initialize some constants
    chainlength <- length(phybreak.object$s$mu)
    obsize <- phybreak.object$p$obs
    samplerange <- (chainlength - samplesize + 1):chainlength
    posteriorsamples <- phybreak.object$s$nodehosts[obsize:(2 * obsize - 1), samplerange]
    
    ### obtaining the result in steps
    
    # matrix containing posterior support for each infector in each tree
    allsupports <- matrix(0, nrow = obsize, ncol = samplesize)
    for (hostID in 1:obsize) {
        for (infector in 0:obsize) {
            allsupports[hostID, posteriorsamples[hostID, ] == infector] <- sum(posteriorsamples[hostID, ] == infector)
        }
    }
    
    # 
    besttree <- which.max(colSums(log(allsupports)))
    if (phylo.class) 
        return(besttree + samplerange[1] - 1)
    res <- matrix(c(posteriorsamples[, besttree], allsupports[, besttree]), ncol = 2)
    
    if (includetimes) {
        timesums <- with(phybreak.object, rowSums(s$nodetimes[obsize:(2 * obsize - 1), samplerange] * (s$nodehosts[obsize:(2 * 
            obsize - 1), samplerange] == res[, 1])))
        timesumsqs <- with(phybreak.object, rowSums((s$nodetimes[obsize:(2 * obsize - 1), samplerange]^2) * (posteriorsamples == 
            res[, 1])))
        res <- cbind(res, timesums/res[, 2])
        res <- cbind(res, sqrt((timesumsqs - res[, 3]^2 * res[, 2])/(res[, 2] - 1)))
        res <- cbind(res, phybreak.object$s$nodetimes[obsize:(2 * obsize - 1), besttree + samplerange[1] - 1])
    }
    
    
    return(res)
}




### function returning the ranknr'th most frequent entry in inf.chain, plus the frequency if support = TRUE called by:
### .transtreecount .inflist (function 'infectorsets')
.postinfector <- function(inf.chain, ranknr = 1, support = FALSE) {
    chainlength <- length(inf.chain)
    ordinf <- rle(sort(inf.chain))
    if (length(ordinf$lengths) < ranknr) {
        if (support) 
            return(c(0, 0)) else return(0)
    }
    ord <- order(ordinf$lengths, decreasing = TRUE)
    if (support) {
        return(c(ordinf[[2]][ord[ranknr]], ordinf[[1]][ord[ranknr]]))
    } else {
        return(ordinf[[2]][ord[ranknr]])
    }
}


### function iteratively removing cycles by edmond's algorithm, as long as there are cycles called by: .transtreeedmonds
### .edmondsiterative calls: .cycleYN .pathtocycle .edmondsiterative
.edmondsiterative <- function(suppmat, samplesize, obs) {
    # get most likely infectors, determine which hosts are in cycles, return ML infectors if there are no cycles
    parentset <- apply(suppmat, 2, which.max)
    cycleIDs <- which(sapply(1:obs, .cycleYN, parset = parentset))
    if (length(cycleIDs) == 0) {
        return(parentset)
    }
    
    # take one cycle, and infector for each member
    cycletonode <- unique(.pathtocycle(parentset, cycleIDs[1]))
    cycleinfector <- parentset[cycletonode]
    
    # determine for each host which cycle member is their most likely infector, and for each host which cycle member would be
    # their infectee
    whichincoming <- cycletonode[apply(suppmat[cycletonode, ], 2, which.max)]
    whichoutgoing <- cycletonode[apply(suppmat[, cycletonode], 1, which.max)]
    
    ## turn the cycle into a single host: in the suppmat, use the position of host #1 for this cycle move all infector support for
    ## any cycle member to host #1 remove all infector support for other hosts let all infectee-candidates be infected by the
    ## cycle (entered into position of host #1)
    suppmat[cycletonode[1], ] <- apply(suppmat[cycletonode, ], 2, max)
    # let all infectee-candidates not be infected by the other positions (removing support)
    suppmat[tail(cycletonode, -1), ] <- -samplesize
    # remove mutual support for cycle members
    suppmat[cycletonode, cycletonode] <- -samplesize
    # adjust the support for infectors of the cycle: the maximum of supports among the cycle members
    suppmat[, cycletonode[1]] <- apply(suppmat[, cycletonode], 1, max) - max(suppmat[, cycletonode])
    # let the other positions all take the index as infector (removing their role in the tree)
    suppmat[, tail(cycletonode, -1)] <- -samplesize
    suppmat[1, tail(cycletonode, -1)] <- 0
    
    # use the adjusted support matrix to get a proper transmission tree
    treeres <- .edmondsiterative(suppmat, samplesize, obs)
    
    # determine which hosts are infected by the cycle, and the infector of the cycle
    incoming <- which(treeres == cycletonode[1])
    outgoing <- treeres[cycletonode[1]]
    
    # let the cycle infectees be infected by their most likely infector
    treeres[incoming] <- whichincoming[incoming]
    # resolve the cycle by choosing the best infectee for the selected infector
    cycleinfector[cycletonode == whichoutgoing[outgoing]] <- outgoing
    treeres[cycletonode] <- cycleinfector
    
    
    return(treeres)
}

### returns TRUE if ID is in a cycle, FALSE if not, based on vector of infectors called by: .edmondsiterative calls:
### .pathtocycle
.cycleYN <- function(parset, ID) {
    .pathtocycle(parset, ID)[1] == ID && ID != 1
}

### return path to the root from IDs[1], iteratively calling this function stop if a cycle is encountered called by:
### .edmondsiterative .cycleYN .pathtocycle calls: .pathtocycle
.pathtocycle <- function(parset, IDs) {
    if (IDs[1] == 1 | length(IDs) > length(unique(IDs))) {
        return(IDs)
    } else {
        return(.pathtocycle(parset, c(parset[IDs[1]], IDs)))
    }
}
