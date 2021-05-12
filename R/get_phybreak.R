#' Extracting from a phybreak object
#' 
#' @param x An object of class \code{phybreak}.
#' 
#' @name get_phybreak
#'   
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First build a phybreak-object.
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' 
#' get_transtree(MCMCstate)
#' get_parameters(MCMCstate)
#' mcmcobject <- get_mcmc(MCMCstate, thin = 2)
#' plot.phylo(get_phylo(MCMCstate))
#' get_data(MCMCstate)
#' get_bottlenecks(MCMCstate)
#' 
#' #function from package phangorn:
#' phangorn::parsimony(get_phylo(MCMCstate), get_data(MCMCstate)$sequences)
#' 
#' tree0 <- get_phylo(MCMCstate)
#' seqdata <- get_data(MCMCstate)$sequences
#' phangorn::pml(tree0, seqdata, 
#'               rate = 0.75*get_parameters(MCMCstate, "mu")) 
#' logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
#'        sampling = FALSE, generation = FALSE) #should give the same result as 'pml'
NULL
#> NULL


### the sequence data in class 'phyDat'
#' @describeIn get_phybreak The data in class \code{phybreakdata}.
#' @export
get_data <- function(x) {
  phybreakdata(sequences = x$d$sequences, sample.times = x$d$sample.times,
               sample.names = x$d$names, host.names = x$d$hostnames)
}

#' @export
get.seqdata <- function(x) {
  .Deprecated("get_data")
  return(x$d$sequences)
}


#' @describeIn get_phybreak A named vector with current parameter values.
#' @param whichpars Which parameters to return. Either a vector with parameter names, or \code{"all"} for all parameters, or 
#'   \code{"posterior"} for parameters for which a posterior is sampled.
#' @export
get_parameters <- function(x, whichpars = "posterior") {
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  if(whichpars == "posterior") {
    whichpars <- c("mu", 
                   if(x$h$est.mG) "gen.mean",
                   if(x$h$est.mS) "sample.mean",
                   if(x$h$est.wh.s) "wh.slope",
                   if(x$h$est.wh.e) "wh.exponent",
                   if(x$h$est.wh.0) "wh.level")
  } else if(whichpars == "all") {
    whichpars <- setdiff(names(x$p), c("obs", "wh.model"))
  } else if(!all(whichpars %in% names(x$p))) {
    warning("some parameters in whichpars don't exist and are left out")
    whichpars <- whichpars[whichpars %in% names(x$p)]
  }
  
  pars <- t(as.data.frame(x$p)[whichpars])[, 1]
  return(pars)
}  

#' @export
get.parameters <- function(x, samplenr = 0, whichpars = "posterior") {
  .Deprecated("get_parameters", msg = "'get.parameters' is deprecated.\nUse 'get_parameters' instead, without argument 'samplenr'.\nSee help(\"Deprecated\")")
  get_parameters(x, whichpars)
}


#' @describeIn get_phybreak A \code{data.frame} with current infectors and infection times.
#' @export
get_transtree <- function(x) {
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  res <- phybreak2trans(x$v, x$d$hostnames, x$d$reference.date)
  res <- with(res, 
              data.frame(infectors = sim.infectors, inf.times = sim.infection.times, row.names = names(sim.infectors)))
  return(res)
}

#' @export
get.tree <- function(x, samplenr = 0) {
  .Deprecated("get_transtree", msg = "'get.tree' is deprecated.\nUse 'get_transtree' instead, without argument 'samplenr'.\nSee help(\"Deprecated\")")
  get_transtree(x)
}


#' @describeIn get_phybreak Returns an object of class \code{\link[ape]{phylo}} ans optionally of class
#'   \code{"simmap"} (package \pkg{phytools}).
#' @param samplenr The posterior tree sample to choose. If \code{samplenr = 0}, the current state is used.
#' @param simmap Whether to include class \code{"simmap"} elements (package \pkg{phytools}), colouring the branches
#'   on the tree to indicate hosts. Is used by \code{\link{plotPhylo}}.
#' @export
get_phylo <- function(x, simmap = FALSE, samplenr = 0) {
  ### tests
  if (simmap & !("phytools" %in% .packages(TRUE))) {
    warning("package 'phytools' is not installed; changed input to simmap = FALSE")
  }
  if (simmap & !("phytools" %in% .packages(FALSE))) {
    warning("package 'phytools' is not attached while simmap = TRUE")
  }
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  if (samplenr > length(x$s$mu)) {
    warning("requested 'samplenr' not available; current state used")
    samplenr <- 0
  }
  
  ### make tree
  if (samplenr == 0) {
    return(phybreak2phylo(x$v, x$d$names, simmap)) 
  } else {
    vars <- with(x, 
                 list(
                   inftimes = s$inftimes[, samplenr],
                   infectors = s$infectors[, samplenr],
                   nodetimes = c(v$nodetimes[v$nodetypes %in% c("s", "x")], s$nodetimes[, samplenr]),
                   nodeparents = s$nodeparents[, samplenr],
                   nodehosts = c(v$nodehosts[v$nodetypes %in% c("s", "x")], s$nodehosts[, samplenr]),
                   nodetypes = v$nodetypes
                 ))
    return(phybreak2phylo(vars, x$d$names, simmap))
  }
}

#' @export
get.phylo <- function(...) {
  .Deprecated("get_phylo")
  get_phylo(...)
}


#' @describeIn get_phybreak Returns an object of class \code{\link[ape]{multiphylo}}.
#'   
#' @export
get_multiPhylo <- function(x, thin = 1, nkeep = Inf) {
  ### tests
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  chainlength <- length(x$s$mu)
  if (nkeep * thin > chainlength & nkeep < Inf) {
    warning("'nkeep * thin' larger than number of available samples. Fewer trees are included.")
  }
  
  ### posterior samples to be included
  tokeep <- seq(thin, chainlength, thin)
  tokeep <- tail(tokeep, nkeep)

  
  ### make trees
  res <- lapply(tokeep, get.phylo, x = x)
  names(res) <- tokeep
  class(res) <- "multiPhylo"
  
  return(res)
}

#' @export
get.multiPhylo <- function(...) {
  .Deprecated("get_multiPhylo")
  get_multiPhylo(...)
}





### an mcmc-object (package 'coda')
#' @describeIn get_phybreak An object of class \code{"mcmc"} (package \pkg{coda}), with sampled parameters, 
#'   infection times, and infectors.
#'   
#' @param thin Thinning interval.
#' @param nkeep Number of samples to keep, counting from tail of the chain.
#' @export
get_mcmc <- function(x, thin = 1, nkeep = Inf) {
  ### tests
  if(!("coda" %in% .packages(TRUE))) {
    stop("package 'coda' should be installed for this function")
  }
  if(!("coda" %in% .packages(FALSE))) {
    warning("package 'coda' is not attached")
  }
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  chainlength <- length(x$s$mu)
  if (nkeep * thin > chainlength & nkeep < Inf) {
    warning("'nkeep * thin' larger than number of available samples. Fewer samples are kept.")
  }
  
  ### posterior samples to be included
  tokeep <- seq(thin, chainlength, thin)
  tokeep <- tail(tokeep, nkeep)
  
  ### extracting all variables and parameters, and naming them
  res <- with(x, cbind(t(s$inftimes[, tokeep]), 
                       t(s$infectors[, tokeep])))
  parnames <- with(x,
                   c(paste0("tinf.", d$hostnames[1:p$obs]), paste0("infector.", d$hostnames[1:p$obs])))
  if (x$h$est.wh.h) {
    res <- cbind(x$s$wh.h[tokeep], res)
    parnames <- c("wh.history", parnames)
  }
  if (x$h$est.wh.s) {
    res <- cbind(x$s$wh.s[tokeep], res)
    parnames <- c("wh.slope", parnames)
  }
  if (x$h$est.wh.e) {
    res <- cbind(x$s$wh.e[tokeep], res)
    parnames <- c("wh.exponent", parnames)
  }
  if (x$h$est.wh.0) {
    res <- cbind(x$s$wh.0[tokeep], res)
    parnames <- c("wh.level", parnames)
  }
  if (x$h$est.mS) {
    res <- cbind(x$s$mS[tokeep], res)
    parnames <- c("mS", parnames)
  }
  if (x$h$est.mG) {
    res <- cbind(x$s$mG[tokeep], res)
    parnames <- c("mG", parnames)
  }
  if (x$h$est.tG) {
    res <- cbind(x$s$tG[tokeep], res)
    parnames <- c("tG", parnames)
  }
  if (x$h$est.tS) {
    res <- cbind(x$s$tS[tokeep], res)
    parnames <- c("tS", parnames)
  }
  res <- cbind(x$s$hist_dens[tokeep], res)
  parnames <- c("hist.dens", parnames)
  res <- cbind(x$s$historyinf[tokeep], res)
  parnames <- c("historyinf", parnames)
  
  res <- cbind(x$s$mu[tokeep], x$s$introductions[tokeep], res, x$s$logLik[tokeep])
  parnames <- c("mu", "introductions", parnames, "logLik")
  colnames(res) <- parnames
  
  return(coda::mcmc(res))
}


#' @export
get.mcmc <- function(...) {
  .Deprecated("get_mcmc")
  get_mcmc(...)
}




### an mcmc-object (package 'coda')
#' @describeIn get_phybreak A vector with for each host the current number of extant lineages in its phylogenetic minitree
#'   at the time of infection (the bottleneck size).
#'   
#' @export
get_bottlenecks <- function(x) {
  v <- phybreak2environment(x$v)
  
  res <- 1L + with(v,
                   tabulate(nodehosts[which(nodeparents %in% which(nodetypes %in% "b"))], x$p$obs))
  
  names(res) <- x$d$hostnames[1:x$p$obs]
  
  return(res)
}

