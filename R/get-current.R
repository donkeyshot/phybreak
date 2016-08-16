### get current state of the tree ### This file contains get.XXX functions

#' Accessing a phybreak object
#' 
#' @param phybreak.object An object of class \code{phybreak}.
#' 
#' @name get.phybreak
NULL
#> NULL

### a 2-column matrix with infectors and infection times
#' @describeIn get.phybreak A \code{data.frame} with current infectors and infection times.
#' @export
get.tree <- function(phybreak.object) {
    if (!inherits(phybreak.object, "phybreak")) {
        stop("object must be of class \"phybreak\"")
    }
    infectors <- c("index", phybreak.object$d$names)[1 + tail(phybreak.object$v$nodehosts, phybreak.object$p$obs)]
    inftimes <- tail(phybreak.object$v$nodetimes, phybreak.object$p$obs)
    res <- data.frame(infectors = infectors, inf.times = inftimes, row.names = phybreak.object$d$names)
    return(res)
}

### a vector with named parameter values
#' @describeIn get.phybreak A named vector with all current parameter values.
#' @export
get.parameters <- function(phybreak.object) {
    return(unlist(phybreak.object$p))
}

### an mcmc-object (package 'coda')
#' @describeIn get.phybreak An object of class \code{"mcmc"} (package \pkg{coda}), with sampled parameters, 
#'   infection times, and infectors.
#'   
#' @param thin Thinning interval for making the mcmc-object.
#' @param nkeep Number of samples to keep in the mcmc-object, counting from tail of the chain.
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
    res <- with(phybreak.object, cbind(t(s$nodetimes[p$obs:(2 * p$obs - 1), tokeep]), t(s$nodehosts[p$obs:(2 * p$obs - 1), tokeep])))
    parnames <- c(paste0("tinf.", phybreak.object$d$names), paste0("infector.", phybreak.object$d$names))
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


### the phylogenetic tree in class 'phylo', and possibly also in class 'simmap'. The function also allows to get one of the
### posterior trees (0 = current) called by: plot.phybreak calls: .makephylosimmap.phybreak .makephylo.phybreak

#' @describeIn get.phybreak Returns an object of class \code{"phylo"} (package \pkg{ape}) and optionally of class 
#'   \code{"simmap"} (package \pkg{phytools}).
#'   
#' @param post.tree The posterior tree sample to choose. If \code{post.tree = 0}, the current state is used.
#' @param simmap Whether to include class \code{"simmap"} elements (package \pkg{phytools}), colouring the branches 
#'   on the tree to indicate hosts.
#'   
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' #First build a phybreak-object.
#' simulation <- sim.phybreak(obsize = 5)
#' MCMCstate <- phybreak(data = simulation$sequences, times = simulation$sample.times)
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
#' phangorn::pml(tree0, seqdata, rate = 0.75*get.parameters(MCMCstate)["mu"]) 
#' logLik(MCMCstate, genetic = TRUE, withinhost = FALSE, 
#'        sampling = FALSE, generation = FALSE) #should give the same result as 'pml'
#' @export
get.phylo <- function(phybreak.object, post.tree = 0, simmap = FALSE) {
    ### tests
    if (simmap & !("phytools" %in% .packages(TRUE))) {
      stop("package 'phytools' should be installed for this function if simmap = TRUE")
    }
    if (simmap & !("phytools" %in% .packages(FALSE))) {
      warning("package 'phytools' is not attached while simmap = TRUE")
    }
    if (!inherits(phybreak.object, "phybreak")) {
        stop("object must be of class \"phybreak\"")
    }
    if (post.tree > length(phybreak.object$s$mu)) {
        warning("requested 'post.tree' not available; current state used")
        post.tree <- 0
    }
    
    ### get times, parents, hosts
    if (post.tree == 0) {
        nodetimes <- phybreak.object$v$nodetimes
        nodeparents <- phybreak.object$v$nodeparents
        nodehosts <- phybreak.object$v$nodehosts
    } else {
        nodetimes <- with(phybreak.object, c(v$nodetimes[1:p$obs], s$nodetimes[, post.tree]))
        nodeparents <- phybreak.object$s$nodeparents[, post.tree]
        nodehosts <- with(phybreak.object, c(v$nodehosts[1:p$obs], s$nodehosts[, post.tree]))
    }
    nodenames <- phybreak.object$d$names
    
    
    if (simmap) {
        return(.makephylosimmap.phybreak(nodetimes, nodeparents, nodehosts, nodenames))
    } else {
        return(.makephylo.phybreak(nodetimes, nodeparents, nodenames))
    }
    
}


### the sequence data in class 'phyDat'
#' @describeIn get.phybreak The sequence data in class \code{"phyDat"} (package \pkg{phangorn}).
#' @export
get.seqdata <- function(phybreak.object) {
    return(phybreak.object$d$sequences)
}



