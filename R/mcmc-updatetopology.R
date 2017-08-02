.updatehost.withinhost <- function(hostID) {
  ### create an up-to-date proposal-environment with hostID as focal host
  .prepare.pbe()
  .copy2pbe1("hostID", environment())
  
  ### remove current state
  rewire_stripminitree(.pbe1, .pbe1$hostID)

  ### change to proposed state
  rewire_buildminitree(.pbe1, .pbe1$hostID)
  
  ### calculate proposal ratio
  logproposalratio <- 0

  ### calculate likelihood
  .propose.pbe("withinhost")
  
  ### calculate acceptance probability
  logaccprob <- .pbe1$logLikseq - .pbe0$logLikseq + logproposalratio
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    .accept.pbe("withinhost")
  }
  
}
