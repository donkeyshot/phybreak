### get current state of the tree ### This file contains get.XXX functions

#' Accessing a phybreak object
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' 
#' @name get.phybreak
NULL
#> NULL

### a 2-column matrix with infectors and infection times
#' @describeIn get.phybreak A \code{data.frame} with current (\code{samplenr = 0}) or sampled infectors and infection times.

#' @export
get.tree <- function(phybreak.object, samplenr = 0) {
  if (!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  if (samplenr == 0) {
    vars <- list(
      nodetimes = phybreak.object$v$nodetimes,
      nodeparents = phybreak.object$v$nodeparents,
      nodehosts = phybreak.object$v$nodehosts,
      nodetypes = phybreak.object$v$nodetypes
    )
  } else {
    # samplenr > 0
    
    vars <- list(
      nodetimes = c(phybreak.object$v$nodetimes[phybreak.object$v$nodetypes == "s"], phybreak.object$s$nodetimes[, samplenr]),
      nodeparents = phybreak.object$s$nodeparents[, samplenr],
      nodehosts = c(phybreak.object$v$nodehosts[phybreak.object$v$nodetypes == "s"], phybreak.object$s$nodehosts[, samplenr]),
      nodetypes = phybreak.object$v$nodetypes
    )
  }
  
  res <- phybreak2trans(vars, phybreak.object$d$hostnames, phybreak.object$d$reference.date)
  infectors <- c("index", phybreak.object$d$names)[1 + phybreak.object$v$nodehosts[phybreak.object$v$nodetypes == "t"]]
  inftimes <- phybreak.object$v$nodetimes[phybreak.object$v$nodetypes == "t"]
  res <- with(res, 
              data.frame(infectors = sim.infectors, inf.times = sim.infection.times, row.names = names(sim.infectors)))
  return(res)
}

### a vector with named parameter values
#' @describeIn get.phybreak A named vector with current (\code{samplenr = 0}) or sampled parameter values.

#' @param whichpars Which parameters to return. Either a vector with parameter names, or \code{"all"} for all parameters, or 
#'   \code{"posterior"} for parameters for which a posterior is sampled.
#' @export
get.parameters <- function(phybreak.object, samplenr = 0, 
                           whichpars = "posterior") {
  if (!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  if(whichpars == "posterior") {
    whichpars <- c("mu", 
                   if(phybreak.object$h$est.mS) "mean.sample",
                   if(phybreak.object$h$est.mG) "mean.gen",
                   if(phybreak.object$h$est.wh) "wh.slope")
  } else if(whichpars == "all") {
    whichpars <- names(phybreak.object$p)
  } else if(!all(whichpars %in% names(phybreak.object$p))) {
    warning("some parameters in whichpars don't exist and are left out")
    whichpars <- whichpars[whichpars %in% names(phybreak.object$p)]
  }
  
  if (samplenr == 0) {
    pars <- unlist(phybreak.object$p)[whichpars]
  } else {
    pars <- c()
    for(param in whichpars) {
      pars <- c(pars, phybreak.object$p[[param]])
      names(pars) <- c(head(names(pars), -1), param)
      if(param == "mu") {
        pars[param] <- phybreak.object$s$mu[samplenr]
      }
      if(param == "mean.sample") {
        pars[param] <- phybreak.object$s$mS[samplenr]
      }
      if(param == "mean.gen") {
        pars[param] <- phybreak.object$s$mG[samplenr]
      }
      if(param == "wh.slope") {
        pars[param] <- phybreak.object$s$slope[samplenr]
      }
    }
  }
  
  return(pars)
}

### an mcmc-object (package 'coda')
#' @describeIn get.phybreak An object of class \code{"mcmc"} (package \pkg{coda}), with sampled parameters, 
#'   infection times, and infectors.
#'   
#' @param thin Thinning interval.
#' @param nkeep Number of samples to keep, counting from tail of the chain.
#' @export
get.mcmc <- function(phybreak.object, thin = 1, nkeep = Inf) {
    ### tests
    if(!("coda" %in% .packages(TRUE))) {
      stop("package 'coda' should be installed for this function")
    }
    if(!("coda" %in% .packages(FALSE))) {
      warning("package 'coda' is not attached")
    }
    if (!inherits(phybreak.object, "phybreak")) {
      stop("object must be of class \"phybreak\"")
    }
    chainlength <- length(phybreak.object$s$mu)
    if (nkeep * thin > chainlength & nkeep < Inf) {
        warning("'nkeep * thin' larger than number of available samples. Fewer samples are kept.")
    }
    
    ### posterior samples to be included
    tokeep <- seq(thin, chainlength, thin)
    tokeep <- tail(tokeep, nkeep)
    
    ### extracting all variables and parameters, and naming them
    res <- with(phybreak.object, cbind(t(s$nodetimes[d$nsamples:(d$nsamples + p$obs - 1), tokeep]), 
                                       t(s$nodehosts[d$nsamples:(d$nsamples + p$obs - 1), tokeep])))
    parnames <- with(phybreak.object,
                     c(paste0("tinf.", d$hostnames[1:p$obs]), paste0("infector.", d$hostnames[1:p$obs])))
    if (phybreak.object$h$est.wh) {
        res <- cbind(phybreak.object$s$slope[tokeep], res)
        parnames <- c("slope", parnames)
    }
    if (phybreak.object$h$est.mG) {
        res <- cbind(phybreak.object$s$mG[tokeep], res)
        parnames <- c("mG", parnames)
    }
    if (phybreak.object$h$est.mS) {
        res <- cbind(phybreak.object$s$mS[tokeep], res)
        parnames <- c("mS", parnames)
    }
    res <- cbind(phybreak.object$s$mu[tokeep], res, phybreak.object$s$logLik[tokeep])
    parnames <- c("mu", parnames, "logLik")
    colnames(res) <- parnames
    
    return(coda::mcmc(res))
}



#' @describeIn get.phybreak Returns an object of class \code{\link[ape]{phylo}} ans optionally of class
#'   \code{"simmap"} (package \pkg{phytools}).
#'   
#' @param samplenr The posterior tree sample to choose. If \code{samplenr = 0}, the current state is used.
#' @param simmap Whether to include class \code{"simmap"} elements (package \pkg{phytools}), colouring the branches
#'   on the tree to indicate hosts. Is used by \code{\link{plotPhylo}}.
#'   
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation)
#' MCMCstate <- burnin.phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample.phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' get.tree(MCMCstate)
#' get.parameters(MCMCstate)
#' codaobject <- get.mcmc(MCMCstate, thin = 2)
#' plot.phylo(get.phylo(MCMCstate))
#' get.seqdata(MCMCstate)
#' 
#' #function from package phangorn:
#' phangorn::parsimony(get.phylo(MCMCstate), get.seqdata(MCMCstate))
#' 
#' tree0 <- get.phylo(MCMCstate)
#' seqdata <- get.seqdata(MCMCstate)
#' phangorn::pml(tree0, seqdata, 
#'               rate = 0.75*get.parameters(MCMCstate)["mu"]) 
#' logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
#'        sampling = FALSE, generation = FALSE) 
#'               #should give the same result as 'pml'
#' @export
get.phylo <- function(phybreak.object, samplenr = 0, simmap = FALSE) {
  ### tests
  if (simmap & !("phytools" %in% .packages(TRUE))) {
    warning("package 'phytools' is not installed; changed input to simmap = FALSE")
  }
  if (simmap & !("phytools" %in% .packages(FALSE))) {
    warning("package 'phytools' is not attached while simmap = TRUE")
  }
  if (!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  if (samplenr > length(phybreak.object$s$mu)) {
    warning("requested 'samplenr' not available; current state used")
    samplenr <- 0
  }
  
  ### make tree
  if (samplenr == 0) {
    return(phybreak2phylo(phybreak.object$v, phybreak.object$d$names, simmap)) 
  } else {
    vars <- with(phybreak.object, 
                 list(
                   nodetimes = c(v$nodetimes[v$nodetypes %in% c("s", "x")], s$nodetimes[, samplenr]),
                   nodeparents = s$nodeparents[, samplenr],
                   nodehosts = c(v$nodehosts[v$nodetypes %in% c("s", "x")], s$nodehosts[, samplenr]),
                   nodetypes = v$nodetypes
                 ))
    return(phybreak2phylo(vars, phybreak.object$d$names, simmap))
  }

}



#' @describeIn get.phybreak Returns an object of class \code{\link[ape]{multiphylo}}.
#'   
#' @export
get.multiPhylo <- function(phybreak.object, thin = 1, nkeep = Inf) {
  ### tests
  if (!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  chainlength <- length(phybreak.object$s$mu)
  if (nkeep * thin > chainlength & nkeep < Inf) {
    warning("'nkeep * thin' larger than number of available samples. Fewer trees are included.")
  }
  
  ### posterior samples to be included
  tokeep <- seq(thin, chainlength, thin)
  tokeep <- tail(tokeep, nkeep)

  
  ### make trees
  res <- lapply(tokeep, get.phylo, phybreak.object = phybreak.object)
  names(res) <- tokeep
  class(res) <- "multiPhylo"
  
  return(res)
}



### the sequence data in class 'phyDat'
#' @describeIn get.phybreak The sequence data in class \code{"phyDat"} (package \pkg{phangorn}).
#' @export
get.seqdata <- function(phybreak.object) {
    return(phybreak.object$d$sequences)
}



