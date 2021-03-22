#' Sampled infectors for each host in a phybreak-object.
#' 
#' The function takes a \code{phybreak}-object containing MCMC-samples, and returns for each host 
#'   the posterior infectors, with support per infector. That can be in the form of a list, with infectors ordered
#'   by support, or in the form of a matrix with supports for all host-infector combinations (host in columns, infectors in rows).
#' 
#' @param x An object of class \code{phybreak}.
#' @param which.hosts A vector with hosts (exact names or positions in the dataset), or \code{"all"} for all hosts.
#' @param percentile Return infectors ordered by support, until a cumulative support indicated by \code{percentile},
#'   support measured by proportion (value between 0 and 1).
#' @param minsupport Only return infectors with more support than \code{minsupport}. Values in the range [0,1) are 
#'   interpreted as support in \code{"proportion"} of posterior samples, values >= 1 as \code{"count"} of posterior 
#'   samples.
#' @param samplesize The number of samples to include (taken from the tail of the MCMC-chain).
#' @param infector.name Whether to return the names of the infectors, or their position in the dataset.
#' @param support Whether to return the support (= posterior probability) for each infector as a \code{"proportion"} 
#'   or as a \code{"count"} of posterior trees in which that transmission link or transmission cluster is present.
#' @param output Whether to return a \code{list} of hosts with infectors, or a \code{matrix} with support values 
#'   for all host-infector pairs.
#'   Arguments \code{percentile}, \code{samplesize}, and \code{infector.name} are not used in case of \code{matrix}.
#' @return A named list with \code{data.frame}s, each with a vector of infectors and a vector of supports. Or a matrix with
#'   supports for each host-infector pair.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object containing samples.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' infectorsets(MCMCstate)
#' @export
infectorsets <- function(x, which.hosts = "all", percentile = 0.95, minsupport = 0, samplesize = Inf, infector.name = TRUE, 
    support = c("proportion", "count"), output = c("list", "matrix")) {
  if (!inherits(x, "phybreak")) 
    stop("object should be of class \"phybreak\"")
  
  ### allow for older phybreak versions
  obs <- x$p$obs
  if(is.null(x$d$nsamples)) {
    nsamples <- obs
  } else {
    nsamples <- x$d$nsamples
  }
  if(is.null(x$d$hostnames)) {
    hostnames <- x$d$names
  } else {
    hostnames <- x$d$hostnames[1:obs]
  }
  if(is.null(x$s$infectors)) {
    infectors <- x$s$nodehosts[nsamples:(nsamples + obs), ]
  } else {
    infectors <- x$s$infectors
  }
  
  ### tests
  if (!(is.character(which.hosts) | is.numeric(which.hosts)) | any(is.na(which.hosts))) {
    stop("'which.hosts' should be numeric or \"all\", or should contain exact host names")
  }
  if(is.numeric(which.hosts)) {
    which.hosts <- as.integer(which.hosts)
    which.hosts <- unique(which.hosts)
    which.hosts <- which.hosts[which.hosts >= 1 & which.hosts <= obs]
  } else if (any(which.hosts == "all")) {
    which.hosts <- 1:obs
  } else {
    which.hosts <- match.arg(which.hosts, hostnames, several.ok = TRUE)
    which.hosts <- unique(which.hosts)
    which.hosts <- match(which.hosts, hostnames)
  } 
  if (percentile < 0 | percentile > 1) {
    stop("'percentile' should be given as number between 0 and 1")
  }
  if (minsupport < 0) {
    minsupport <- 0
    warning("'minsupport' is set to 0 (it should be >= 0). Values in the range [0,1] are interpreted as support in \"proportion\" of 
        posterior samples, values > 1 as \"count\" of posterior samples" )
  }
  if (ncol(infectors) == 0) stop("no sampled trees available")
  if (samplesize > ncol(infectors) & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  support <- match.arg(support)
  output <- match.arg(output)

    ### samplesize
    chainlength <- ncol(infectors)
    samplesize <- min(samplesize, chainlength)
    samplerange <- (chainlength - samplesize + 1):chainlength
    if (support == "count") {
      denominator <- 1
    } else {
      denominator <- samplesize
    }
    

    ### obtaining the result in steps
    
    # array with ordered infectors per host, plus support
    
    if(output == "matrix") {
      supportmatrix <- apply(1 + infectors[, samplerange], 
                             1, tabulate, nbin = 1 + obs)
      colnames(supportmatrix) <- hostnames
      rownames(supportmatrix) <- c("index", hostnames)
      return(supportmatrix[, which.hosts]/denominator)
    }
    
    inffreqs <- .infarray(infectors[, samplerange])
    
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
            res <- c(res, list(data.frame(infector = c("index", hostnames)[1 + inffreqs[i, 1, ][includemat[i, ]]], 
                support = inffreqs[i, 2, ][includemat[i, ]]/denominator, row.names = NULL)))
        }
    } else {
        res <- list()
        for (i in which.hosts) {
            res <- c(res, list(data.frame(infector = inffreqs[i, 1, ][includemat[i, ]], support = inffreqs[i, 2, ][includemat[i, 
                ]]/denominator, row.names = NULL)))
        }
    }
    names(res) <- hostnames[which.hosts]
    
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


