### run mcmc chain and return the updated phylo-object ###

### phybreak functions called ###
# .build.phybreakenv
# .updatehost
# .updatehost.keepphylo
# .update.mG
# .update.mS
# .update.wh
# .update.mu
# .destroy.phybreakenv


burnin.phybreak <- function(phybreak.object, ncycles, keepphylo = 0.2) {
  .build.phybreakenv(phybreak.object)
  
  
  for(rep in 1:ncycles) {
    for(i in sample(phybreak.object$p$obs)) {
      if(runif(1) < 1 - keepphylo) .updatehost(i) else .updatehost.keepphylo(i)
    }
    if(phybreak.object$h$est.mG) .update.mG()
    if(phybreak.object$h$est.mS) .update.mS()
    if(phybreak.object$h$est.wh) .update.wh()
    .update.mu()
  }
  
  res <- .destroy.phybreakenv(phybreak.object$s)
  

  return(res)
}
