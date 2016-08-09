### get current state of the tree ###
### This file contains get.XXX functions 

### a 1-column matrix with infectors
get.infectors <- function(phybreak.object) {
  if(!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  infectors <- c("index", phybreak.object$d$names)[
    1 + tail(phybreak.object$v$nodehosts,phybreak.object$p$obs)]
  res <- matrix(infectors, ncol = 1,
                dimnames = list(phybreak.object$d$names, "infector"))
  return(res)
}

### a vector with named parameter values
get.parameters <- function(phybreak.object) {
  return(unlist(phybreak.object$p))
}

### an mcmc-object (package 'coda')
get.mcmc <- function(phybreak.object, thin = 1, samplesize = Inf) {
  ### tests
  if(!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  chainlength <- length(phybreak.object$s$mu)
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  
  ### posterior samples to be included
  tokeep <- seq(thin, length(phybreak.object$s$logLik), thin)
  tokeep <- tail(tokeep, samplesize)
  
  ### extracting all variables and parameters, and naming them
  res <- with(phybreak.object,
              cbind(t(s$nodetimes[p$obs:(2*p$obs-1), tokeep]),
                    t(s$nodehosts[p$obs:(2*p$obs-1), tokeep])))
  parnames <- c(paste0("tinf.",phybreak.object$d$names),
                paste0("infector.",phybreak.object$d$names))
  if(phybreak.object$h$est.wh) {
    res <- cbind(phybreak.object$s$slope[tokeep], res)
    parnames <- c("slope",parnames)
  }
  if(phybreak.object$h$est.mG) {
    res <- cbind(phybreak.object$s$mG[tokeep], res)
    parnames <- c("mG",parnames)
  }
  if(phybreak.object$h$est.mS) {
    res <- cbind(phybreak.object$s$mS[tokeep], res)
    parnames <- c("mS",parnames)
  }
  res <- cbind(phybreak.object$s$mu[tokeep], res,
               phybreak.object$s$logLik[tokeep])
  parnames <- c("mu",parnames,"logLik")
  colnames(res) <- parnames
  
  return(mcmc(res))
}


### the phylogenetic tree in class 'phylo', and possibly also
### in class 'simmap'. The function also allows to get one of
### the posterior trees (0 = current) 
### called by:
# plot.phybreak
### calls:
# .makephylosimmap.phybreak
# .makephylo.phybreak
get.phylo <- function(phybreak.object, post.tree = 0, simmap = FALSE) {
  ### tests
  if(!inherits(phybreak.object, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  if(post.tree > length(phybreak.object$s$mu)) {
    warning("requested 'post.tree' not available; current state used")
    post.tree <- 0
  }
  
  ### get times, parents, hosts
  if(post.tree == 0) {
    nodetimes <- phybreak.object$v$nodetimes
    nodeparents <- phybreak.object$v$nodeparents
    nodehosts <- phybreak.object$v$nodehosts
  } else {
    nodetimes <- with(phybreak.object,
                      c(v$nodetimes[1:p$obs],
                        s$nodetimes[,post.tree]))
    nodeparents <- phybreak.object$s$nodeparents[,post.tree]
    nodehosts <- with(phybreak.object,
                      c(v$nodehosts[1:p$obs],
                        s$nodehosts[,post.tree]))
  }
  nodenames <- phybreak.object$d$names
  
  
  if(simmap) {
    return(.makephylosimmap.phybreak(nodetimes, nodeparents, nodehosts, nodenames))
  } else {
    return(.makephylo.phybreak(nodetimes, nodeparents, nodenames))
  }
  
}
  

### the sequence data in class 'phyDat' 
get.phyDat <- function(phybreak.object) {
  ###sequences
  dnadata <- with(phybreak.object, {
    t(matrix(rep(t(d$SNP),rep(d$SNPfr,p$obs)),ncol=p$obs))
  })
  rownames(dnadata) <- phybreak.object$d$names
  
  return(phyDat(dnadata))
}



