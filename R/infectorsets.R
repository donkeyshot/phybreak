### sets of possible infectors per host, from phybreak-object with samples ###

### function returns a list with for each host their most likely infectors (ordered)
### with support. 
### calls:
# .infarray
infectorsets <- function(phybreak.object, which.hosts = "all", percentile = 0.95, minsupport = 0, 
                         samplesize = Inf, infector.name = TRUE, support = "proportion"
                         ) {
  ### initialize some constants
  chainlength <- length(phybreak.object$s$mu)
  obs <- phybreak.object$p$obs
  samplesize <- min(samplesize, chainlength)
  samplerange <- (chainlength-samplesize+1):chainlength
  if(which.hosts == "all") {which.hosts <- 1:obs}
  if(support == "count") {denominator <- 1} else {denominator <- samplesize}

  ### tests
  if(chainlength == 0) stop("no sampled trees available")
  if(samplesize > chainlength & samplesize < Inf) {
    warning("desired 'samplesize' larger than number of available samples")
  }
  if(support != "proportion" && support != "count") {
    warning("support is given as proportion")
  }
  if(!is.numeric(which.hosts) || !is.integer(which.hosts) || any(is.na(which.hosts))) {
    which.hosts <- match(which.hosts, phybreak.object$d$names)
    if(!is.numeric(which.hosts) || !is.integer(which.hosts) || any(is.na(which.hosts))) {
      stop("which.hosts should be numeric or \"all\", or should contain exact host names")
    }
  }

  ### obtaining the result in steps

  # array with ordered infectors per host, plus support  
  inffreqs <- .infarray(phybreak.object$s$nodehosts[obs:(2*obs-1),samplerange])
  
  # take out those with too little support
  includemat <- inffreqs[,2,] > minsupport * samplesize
  
  # only include up to a total cumulative support
  cumfreqincl <- t(rbind(rep(TRUE, obs),
                         apply(inffreqs[,2,],1,cumsum) < percentile * samplesize)[-(obs+1),])
  includemat <- includemat & cumfreqincl
  
  ### construct the output
  if(infector.name) {
    res <- list()
    for(i in which.hosts) {
      res <- c(res,list(data.frame(
        infector = c("index", phybreak.object$d$names)[1 + inffreqs[i,1,][includemat[i,]]],
        support = inffreqs[i,2,][includemat[i,]]/denominator)))
    }
  } else {
    res <- list()
    for(i in which.hosts) {
      res <- c(res,list(data.frame(
        infector = inffreqs[i,1,][includemat[i,]],
        support = inffreqs[i,2,][includemat[i,]]/denominator)))
    }
  }
  names(res) <- phybreak.object$d$names
  
  ### return the result
  return(res)
}


### this function returns ordered most likely infectors for multiple hosts, with posterior support
### called from:
# infectorsets
### calls:
# .inflist
.infarray <- function(inf.matrix) {
  res <- t(apply(inf.matrix,MARGIN = 1, FUN = .inflist, nhosts = nrow(inf.matrix)))
  dim(res) <- c(nrow(inf.matrix), 2, nrow(inf.matrix))
  return(res)
}

### this function returns ordered most likely infectors for a single host, with posterior support
### called from:
# .infarray
### calls:
# .postinfector (file 'transtree')
.inflist <- function(inf.chain, nhosts) {
  sapply(1:nhosts, .postinfector, inf.chain=inf.chain, support = TRUE)
}

