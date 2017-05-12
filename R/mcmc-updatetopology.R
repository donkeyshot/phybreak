.updatehost.topology <- function(hostID) {
  ### create an up-to-date proposal-environment with hostID as focal host
  .prepare.pbe()
  .copy2pbe1("hostID", environment())
  
  ### making variables and parameters available within the function
  le <- environment()
  d <- .pbe0$d
  p <- .pbe1$p
  v <- .pbe1$v
  
  # phylotree in hostID
  v$nodeparents[v$nodehosts == hostID] <- 
    .sampletopology(which(v$nodehosts == hostID), 
                    v$nodetimes[v$nodehosts == hostID], 
                    v$nodetypes[v$nodehosts == hostID], 
                    hostID + 2 * d$nsamples - 1, p$wh.model)
  
  ### update proposal environment
  .copy2pbe1("v", le)
  
  ### calculate proposal ratio
  logproposalratio <- 0

  ### calculate likelihood
  .propose.pbe("topology")
  
  ### calculate acceptance probability
  logaccprob <- .pbe1$logLikseq - .pbe0$logLikseq + logproposalratio
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    .accept.pbe("topology")
  }
  
}
