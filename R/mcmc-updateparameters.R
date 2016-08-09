### all helper functions exclusively involved in updating the parameters ###


### updating mu by proposing within .phybreakenv.prop,
### and accepting or rejecting
### called from:
# burnin.phybreak
# sample.phybreak
### calling:
# .prepare.phybreakenv
# .propose.phybreakenv
# .accept.phybreakenv
.update.mu <- function() {
  ### create an up-to-date proposal-environment
  .prepare.phybreakenv()
  
  ### change to proposal state
  evalq(
    p$mu <- exp(log(p$mu) +
                rnorm(1,0,h$si.mu)),
    .phybreakenv.prop
  )
  
  ### calculate likelihood
  .propose.phybreakenv("mu")
  
  ### calculate acceptance probability
  evalq(
    logaccprob <- logLik - .phybreakenv$logLik,
    .phybreakenv.prop)
  
  ### accept or reject
  if(runif(1) < exp(.phybreakenv.prop$logaccprob)) {
    .accept.phybreakenv()
  }
}


### updating mS by sampling from the posterior,
### conditional on the current tree
### called from:
# burnin.phybreak
# sample.phybreak
### calling:
# .prepare.phybreakenv
# .propose.phybreakenv
# .accept.phybreakenv
.update.mS <- function() {
  ### create an up-to-date proposal-environment
  .prepare.phybreakenv()
  
  ### change to proposal state
  evalq({
    sumst <- sum(v$nodetimes[v$nodetypes == "s"] -
                   v$nodetimes[v$nodetypes == "t"])
    p$mean.sample <- p$shape.sample / 
      rgamma(1,
             shape = p$shape.sample * p$obs + 2 + (h$mS.av/h$mS.sd)^2,
             rate = sumst + (h$mS.av/p$shape.sample)*(1 + (h$mS.av/h$mS.sd)^2))
  },
  .phybreakenv.prop)
  
  ### calculate likelihood
  .propose.phybreakenv("mS")

  ### accept
  .accept.phybreakenv()
}


### updating mG by sampling from the posterior,
### conditional on the current tree
### called from:
# burnin.phybreak
# sample.phybreak
### calling:
# .prepare.phybreakenv
# .propose.phybreakenv
# .accept.phybreakenv
.update.mG <- function() {
  ### create an up-to-date proposal-environment
  .prepare.phybreakenv()
  
  ### change to proposal state
  evalq({
    sumgt <- sum(v$nodetimes[v$nodetypes == "t" & v$nodehosts != 0] -
                   v$nodetimes[v$nodetypes == "t"][v$nodehosts[v$nodetypes == "t"]])
    p$mean.gen <- p$shape.gen / 
      rgamma(1,
             shape = p$shape.gen * (p$obs-1) + 2 + (h$mG.av/h$mG.sd)^2,
             rate = sumgt + (h$mG.av/p$shape.gen)*(1 + (h$mG.av/h$mG.sd)^2))
  },
  .phybreakenv.prop)
  
  ### calculate likelihood
  .propose.phybreakenv("mG")
  
  ### accept
  .accept.phybreakenv()
}


### updating slope by proposing within .phybreakenv.prop,
### and accepting or rejecting
### called from:
# burnin.phybreak
# sample.phybreak
### calling:
# .prepare.phybreakenv
# .propose.phybreakenv
# .accept.phybreakenv
.update.wh <- function() {
  ### create an up-to-date proposal-environment
  .prepare.phybreakenv()

  evalq({
    ### change to proposal state
    p$wh.slope <- exp(log(p$wh.slope) +
                          rnorm(1,0,h$si.wh))

    ### calculate proposalratio
    logproposalratio <- 
      log(p$wh.slope) - log(.phybreakenv$p$wh.slope)
    },
  .phybreakenv.prop)
  
  ### calculate likelihood
  .propose.phybreakenv("slope")
  
  ### calculate acceptance probability
  evalq(
    logaccprob <- logLikcoal - .phybreakenv$logLikcoal + logproposalratio +
      dgamma(p$wh.slope, shape = h$wh.sh, 
             scale = h$wh.av/h$wh.sh, log = TRUE) - 
      dgamma(.phybreakenv$p$wh.slope, shape = h$wh.sh, 
             scale = h$wh.av/h$wh.sh, log = TRUE),
    .phybreakenv.prop)
  
  ### accept or reject
  if(runif(1) < exp(.phybreakenv.prop$logaccprob)) {
    .accept.phybreakenv()
  }
}
