### sets of possible infectors per host, from phybreak-object with samples ###

### function returns a list with for each host their most likely infectors (ordered) with support.  calls: .infarray

#' Sampled infectors for each host in a phybreak-object.
#' 
#' The function takes a \code{phybreak}-object containing MCMC-samples, and returns for each host 
#'   a table with posterior infectors, with support per infector.
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' @param which.hosts A vector with hosts (positions in the dataset), or \code{"all"} for all hosts.
#' @param percentile Return infectors ordered by support, until a cumulative support indicated by \code{percentile},
#'   support measured by proportion (value between 0 and 1).
#' @param minsupport Only return infectors with more support than \code{minsupport}. Values in the range [0,1] are 
#'   interpreted as support in \code{"proportion"} of posterior samples, values > 1 as code{"count"} of posterior 
#'   samples.
#' @param samplesize The number of samples to include (taken from the tail of the MCMC-chain).
#' @param infector.name Whether to return the names of the infectors, or their position in the dataset.
#' @param support Whether to return the support (= posterior probability) for each infector as a \code{"proportion"} 
#'   or as a \code{"count"} of posterior trees in which that transmission link or transmission cluster is present.
#' @return A named list with \code{data.frame}s, each with a vector of infectors and a vector of supports.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' infectorsets(MCMCstate)
#' @export
infectorsets <- function(phybreak.object, which.hosts = "all", percentile = 0.95, minsupport = 0, samplesize = Inf, infector.name = TRUE, 
    support = c("proportion", "count")) {
    ### initialize some constants
    chainlength <- length(phybreak.object$s$mu)
    obs <- phybreak.object$p$obs
    samplesize <- min(samplesize, chainlength)
    samplerange <- (chainlength - samplesize + 1):chainlength
    if (support[1] == "count") {
      denominator <- 1
    } else {
      denominator <- samplesize
    }
    
    ### tests
    if (!(is.character(which.hosts) | is.numeric(which.hosts)) | any(is.na(which.hosts))) {
      stop("'which.hosts' should be numeric or \"all\", or should contain exact host names")
    }
    which.hosts <- unique(which.hosts)
    if(is.numeric(which.hosts)) {
      if(min(which.hosts) < 1 | max(which.hosts) > phybreak.object$p$obs) {
        stop("'which.hosts' contains out-of-range numbers, should be between 1 and outbreak size")
      }
    } else if (any(which.hosts == "all")) {
      which.hosts <- 1:obs
    } else if (!all(which.hosts %in% phybreak.object$d$names)) {
      stop("'which.hosts' contains non-existing host names'")
    } else which.hosts <- match(which.hosts, phybreak.object$d$names)
    if (percentile < 0 | percentile > 1) {
      stop("'percentile' should be given as number between 0 and 1")
    }
    if (minsupport < 0) {
      stop("'minsupport' should be >= 0. Values in the range [0,1] are interpreted as support in \"proportion\" of 
        posterior samples, values > 1 as \"count\" of posterior samples" )
    }
    if (chainlength == 0) stop("no sampled trees available")
    if (samplesize > chainlength & samplesize < Inf) {
        warning("desired 'samplesize' larger than number of available samples")
    }
    if (support[1] != "proportion" && support[1] != "count") {
        warning("support is given as proportion")
    }

    ### obtaining the result in steps
    
    # array with ordered infectors per host, plus support
    inffreqs <- .infarray(phybreak.object$s$nodehosts[obs:(2 * obs - 1), samplerange])
    
    # take out those with too little support
    if(minsupport < 1) {
      includemat <- inffreqs[, 2, ] > minsupport * samplesize
    } else {
      includemat <- inffreqs[, 2, ] > minsupport
    }

    # only include up to a total cumulative support
    cumfreqincl <- t(rbind(rep(TRUE, obs), apply(inffreqs[, 2, ], 1, cumsum) < percentile * samplesize)[-(obs + 1), ])
    includemat <- includemat & cumfreqincl
    
    ### construct the output
    if (infector.name) {
        res <- list()
        for (i in which.hosts) {
            res <- c(res, list(data.frame(infector = c("index", phybreak.object$d$names)[1 + inffreqs[i, 1, ][includemat[i, ]]], 
                support = inffreqs[i, 2, ][includemat[i, ]]/denominator)))
        }
    } else {
        res <- list()
        for (i in which.hosts) {
            res <- c(res, list(data.frame(infector = inffreqs[i, 1, ][includemat[i, ]], support = inffreqs[i, 2, ][includemat[i, 
                ]]/denominator)))
        }
    }
    names(res) <- phybreak.object$d$names[which.hosts]
    
    ### return the result
    return(res)
}


### this function returns ordered most likely infectors for multiple hosts, with posterior support called from: infectorsets
### calls: .inflist
.infarray <- function(inf.matrix) {
    res <- t(apply(inf.matrix, MARGIN = 1, FUN = .inflist, nhosts = nrow(inf.matrix)))
    dim(res) <- c(nrow(inf.matrix), 2, nrow(inf.matrix))
    return(res)
}

### this function returns ordered most likely infectors for a single host, with posterior support called from: .infarray calls:
### .postinfector (file 'transtree')
.inflist <- function(inf.chain, nhosts) {
    sapply(1:nhosts, .postinfector, inf.chain = inf.chain, support = TRUE)
}

