### all helper functions exclusively involved in updating the parameters ###


### updating mu by proposing within .pbe1, and accepting or rejecting called from: burnin.phybreak sample.phybreak
### calling: .prepare.pbe .propose.pbe .accept.pbe
.update.mu <- function() {
    ### create an up-to-date proposal-environment
    .prepare.pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- .pbe0$h
    p <- .pbe1$p

    ### change to proposal state
    p$mu <- exp(log(p$mu) + rnorm(1, 0, h$si.mu))
    
    ### update proposal environment
    .copy2pbe1("p", le)
    
    ### calculate likelihood
    .propose.pbe("mu")
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikseq - .pbe0$logLikseq
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
        .accept.pbe("mu")
    }
}


### updating mS by sampling from the posterior, conditional on the current tree called from: burnin.phybreak sample.phybreak
### calling: .prepare.pbe .propose.pbe .accept.pbe
.update.mS <- function() {
    ### create an up-to-date proposal-environment
    .prepare.pbe()

    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
    
    ### change to proposal state
    sumst <- sum(v$nodetimes[v$nodetypes == "s"] - v$inftimes)
    p$mean.sample <- p$shape.sample/rgamma(1, shape = p$shape.sample * p$obs + 2 + (h$mS.av/h$mS.sd)^2, rate = sumst + (h$mS.av/p$shape.sample) * 
         (1 + (h$mS.av/h$mS.sd)^2))

    ### update proposal environment
    .copy2pbe1("p", le)
    
    ### calculate likelihood
    .propose.pbe("mS")
    
    ### accept
    .accept.pbe("mS")
}


### updating mG by sampling from the posterior, conditional on the current tree called from: burnin.phybreak sample.phybreak
### calling: .prepare.pbe .propose.pbe .accept.pbe
.update.mG <- function() {
    ### create an up-to-date proposal-environment
    .prepare.pbe()

    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
  
    ### change to proposal state
    sumgt <- sum(v$inftimes[v$infectors != 0] - 
                   v$inftimes[v$infectors])
    p$mean.gen <- p$shape.gen/rgamma(1, 
                                     shape = p$shape.gen * (p$obs - 1) + 2 + (h$mG.av/h$mG.sd)^2, 
                                     rate = sumgt + (h$mG.av/p$shape.gen) * (1 + (h$mG.av/h$mG.sd)^2))
    
    ### update proposal environment
    .copy2pbe1("p", le)
    
    ### calculate likelihood
    .propose.pbe("mG")
    
    ### accept
    .accept.pbe("mG")
}


### updating slope by proposing within .pbe1, and accepting or rejecting called from: burnin.phybreak
### sample.phybreak calling: .prepare.pbe .propose.pbe .accept.pbe
.update.wh.slope <- function() {
    ### create an up-to-date proposal-environment
    .prepare.pbe()

    ### making variables and parameters available within the function
    le <- environment()
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v

    ### change to proposal state
    p$wh.slope <- exp(log(p$wh.slope) + rnorm(1, 0, h$si.wh))

    ### update proposal environment
    .copy2pbe1("p", le)

    ### calculate proposalratio
    logproposalratio <- log(p$wh.slope) - log(.pbe0$p$wh.slope)

    ### calculate likelihood
    .propose.pbe("slope")
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikcoal - .pbe0$logLikcoal + logproposalratio + 
      dgamma(.pbe1$p$wh.slope, shape = h$wh.s.sh, scale = h$wh.s.av/h$wh.s.sh, log = TRUE) - 
      dgamma(.pbe0$p$wh.slope, shape = h$wh.s.sh, scale = h$wh.s.av/h$wh.s.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
        .accept.pbe("slope")
    }
}

.update.wh.exponent <- function() {
  ### create an up-to-date proposal-environment
  .prepare.pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- .pbe0$h
  p <- .pbe1$p
  v <- .pbe1$v
  
  ### change to proposal state
  p$wh.exponent <- exp(log(p$wh.exponent) + rnorm(1, 0, h$si.wh))
  
  ### update proposal environment
  .copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$wh.exponent) - log(.pbe0$p$wh.exponent)
  
  ### calculate likelihood
  .propose.pbe("exponent")
  
  ### calculate acceptance probability
  logaccprob <- .pbe1$logLikcoal - .pbe0$logLikcoal + logproposalratio + 
    dgamma(.pbe1$p$wh.exponent, shape = h$wh.e.sh, scale = h$wh.e.av/h$wh.e.sh, log = TRUE) - 
    dgamma(.pbe0$p$wh.exponent, shape = h$wh.e.sh, scale = h$wh.e.av/h$wh.e.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    .accept.pbe("exponent")
  }
}

.update.wh.level <- function() {
  ### create an up-to-date proposal-environment
  .prepare.pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- .pbe0$h
  p <- .pbe1$p
  v <- .pbe1$v
  
  ### change to proposal state
  p$wh.level <- exp(log(p$wh.level) + rnorm(1, 0, h$si.wh))
  
  ### update proposal environment
  .copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$wh.level) - log(.pbe0$p$wh.level)
  
  ### calculate likelihood
  .propose.pbe("level")
  
  ### calculate acceptance probability
  logaccprob <- .pbe1$logLikcoal - .pbe0$logLikcoal + logproposalratio + 
    dgamma(.pbe1$p$wh.level, shape = h$wh.0.sh, scale = h$wh.0.av/h$wh.0.sh, log = TRUE) - 
    dgamma(.pbe0$p$wh.level, shape = h$wh.0.sh, scale = h$wh.0.av/h$wh.0.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    .accept.pbe("level")
  }
}
