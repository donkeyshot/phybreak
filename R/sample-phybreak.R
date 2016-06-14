### run mcmc chain, take samples from the chain, and return the updated phylo-object ###

### phybreak functions called ###
# .build.phybreakenv
# .updatehost
# .updatehost.keepphylo
# .update.mG
# .update.mS
# .update.wh
# .update.mu
# .destroy.phybreakenv


sample.phybreak <- function(phybreak.object, nsample, thin, keepphylo = 0.2) {
  
  ### create room in s to add the new posterior samples
  s.post <- list(
    nodetimes = with(phybreak.object,cbind(s$nodetimes, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodehosts = with(phybreak.object,cbind(s$nodehosts, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodeparents = with(phybreak.object,cbind(s$nodeparents, matrix(NA, nrow=3*p$obs - 1, ncol=nsample))),
    mu = c(phybreak.object$s$mu, rep(NA, nsample)),
    mG = c(phybreak.object$s$mG, rep(NA, nsample)),
    mS = c(phybreak.object$s$mS, rep(NA, nsample)),
    slope = c(phybreak.object$s$slope, rep(NA, nsample)),
    logLik = c(phybreak.object$s$logLik, rep(NA, nsample))
  )
  
  .build.phybreakenv(phybreak.object)
  
  for(sa in tail(1:length(s.post$mu), nsample)) {
    for(rep in 1:thin) {
      for(i in sample(phybreak.object$p$obs)) {
        if(runif(1) < 1 - keepphylo) .updatehost(i) else .updatehost.keepphylo(i)
      }
      if(phybreak.object$h$est.mG) .update.mG()
      if(phybreak.object$h$est.mS) .update.mS()
      if(phybreak.object$h$est.wh) .update.wh()
      .update.mu()
    }
    s.post$nodetimes[, sa] <- tail(.phybreakenv$v$nodetimes, -phybreak.object$p$obs)
    s.post$nodehosts[, sa] <- tail(.phybreakenv$v$nodehosts, -phybreak.object$p$obs)
    s.post$nodeparents[, sa] <- .phybreakenv$v$nodeparents
    s.post$mu[sa] <- .phybreakenv$p$mu
    s.post$mG[sa] <- .phybreakenv$p$mean.gen
    s.post$mS[sa] <- .phybreakenv$p$mean.sample
    s.post$slope[sa] <- .phybreakenv$p$wh.slope
    s.post$logLik[sa] <- .phybreakenv$logLik + .phybreakenv$logLiksam + 
      + .phybreakenv$logLikgen + .phybreakenv$logLikcoal
  }
  
  res <- .destroy.phybreakenv(s.post)
  
  
  return(res)
  
}
