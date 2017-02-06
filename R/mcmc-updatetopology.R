.updatehost.topology <- function(hostID) {
  ### create an up-to-date proposal-environment with hostID as focal host
  .prepare.pbe()
  .copy2pbe1("hostID", environment())
  
  ### making variables and parameters available within the function
  le <- environment()
  p <- .pbe1$p
  v <- .pbe1$v
  Ngenes <- dim(v$nodetimes)[1]
  
  # phylotree in hostID
  if (v$reassortment[hostID] == 0) { 
    v$nodeparents[, v$nodehosts == hostID] <- 
    .sampletopology(which(v$nodehosts == hostID), 
                    v$nodetimes[1, v$nodehosts == hostID], 
                    v$nodetypes[v$nodehosts == hostID], 
                    hostID + 2 * p$obs - 1, p$wh.model)
  } else {
    v$nodeparents[ ,v$nodehosts == hostID] <- 
      t(sapply(1:Ngenes,
               function(x) .sampletopology(which(v$nodehosts == hostID),
                                           v$nodetimes[x , v$nodehosts == hostID],
                                           v$nodetypes[v$nodehosts == hostID],
                                           hostID + 2 * p$obs - 1, p$wh.model)))
  }
  
  ### update proposal environment
  .copy2pbe1("v", le)
  
  ### calculate proposal ratio
  logproposalratio <- 0

  ### calculate likelihood
  .propose.pbe("topology", hosts = hostID)
  
  ### calculate acceptance probability
  logaccprob <- sum(.pbe1$logLikseq) - sum(.pbe0$logLikseq) + logproposalratio
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    .accept.pbe("topology")
  }
  
}
