### functions exclusively involved in updating the parameters ###


update_mu <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p

    ### change to proposal state
    p$mu <- exp(log(p$mu) + rnorm(1, 0, h$si.mu))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate likelihood
    propose_pbe("mu")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq - pbe0$logLikseq +
      dnorm(log10(pbe1$p$mu), mean = h$mu.av, sd = h$mu.sd, log = TRUE) - 
      dnorm(log10(pbe0$p$mu), mean = h$mu.av, sd = h$mu.sd, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
        accept_pbe("mu")
    }
}


update_mS <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()

    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    sumst <- sum(v$nodetimes[v$nodetypes == "s"] - v$inftimes[v$nodehosts[v$nodetypes == "s"]])
    p$sample.mean <- p$sample.shape/rgamma(1, shape = p$sample.shape * p$obs + 2 + (h$mS.av/h$mS.sd)^2, rate = sumst + (h$mS.av/p$sample.shape) * 
         (1 + (h$mS.av/h$mS.sd)^2))

    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate likelihood
    propose_pbe("mS")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("mS")
    }
}

update_mG <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()

    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
  
    ### change to proposal state
    sumgt <- sum(v$inftimes[v$infectors != 0] - 
                   v$inftimes[v$infectors])
    nrof_indices <- sum(v$infectors == 0)
    p$gen.mean <- p$gen.shape/rgamma(1, 
                                     shape = p$gen.shape * (p$obs - nrof_indices) + 2 + (h$mG.av/h$mG.sd)^2, 
                                     rate = sumgt + (h$mG.av/p$gen.shape) * (1 + (h$mG.av/h$mG.sd)^2))
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate likelihood
    propose_pbe("mG")
    
    ### accept
    accept_pbe("mG")
}

update_tG <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  d <- pbe0$d
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  ### change to proposal state
  sumgt <- sum(v$inftimes[v$infectors > 1] -
                 v$inftimes[v$infectors[v$infectors>1]])
  # thalf <- p$gen.shape/rgamma(1,
  #                                  shape = p$gen.shape * (p$obs - 1) + 2 + (h$mG.av/h$mG.sd)^2,
  #                                  rate = sumgt + (h$mG.av/p$gen.shape) * (1 + (h$mG.av/h$mG.sd)^2))
  
  #p$trans.growth <- max(0,rnorm(1, mean = log(sumgt)/1, sd = 0.1))#log((1-p$trans.init)/p$trans.init)/rexp(1, rate = 15/sumgt)
  p$trans.growth <- exp(log(p$trans.growth) + rnorm(1, 0, 0.1))
  ### update proposal environment
  copy2pbe1("p", le)
  ### calculate likelihood
  propose_pbe("mG")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikgen - pbe0$logLikgen
  
  ### accept
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("mG")
  }
}

update_tS <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  p <- pbe1$p
  v <- pbe1$v
  
  p$trans.sample <- exp(log(p$trans.sample) + rnorm(1, 0, 0.1))
  # 
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate likelihood
  propose_pbe("mG")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikgen - pbe0$logLikgen
  
  ### accept
  if (runif(1) < exp(logaccprob) & p$trans.sample <= 1) {
    accept_pbe("mG")
  }
}

update_ir <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  p$intro.rate <- exp(log(p$intro.rate) + rnorm(1, 0, h$si.ir))
  # 
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$intro.rate) - log(pbe0$p$intro.rate)
  
  ### calculate likelihood
  propose_pbe("ir")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikgen - pbe0$logLikgen
  
  ### accept
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("ir")
  }
}

update_wh_history <- function(){
    ### create an up-to-date proposal-environment
    prepare_pbe()
    
    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    
    ### change to proposal state
    p$wh.history <- exp(log(p$wh.history) + rnorm(1, 0, h$si.wh))
    #if (p$wh.history > 1) return()
    
    ### update proposal environment
    copy2pbe1("p", le)
    
    ### calculate proposalratio
    logproposalratio <- log(p$wh.history) - log(pbe0$p$wh.history)
    
    ### calculate likelihood
    propose_pbe("wh.history")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal + logproposalratio + 
      dgamma(pbe1$p$wh.history, shape = h$wh.h.sh, scale = h$wh.h.av/h$wh.h.sh, log = TRUE) - 
      dgamma(pbe0$p$wh.history, shape = h$wh.h.sh, scale = h$wh.h.av/h$wh.h.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("wh.history")
    }
}

update_wh_slope <- function() {
    ### create an up-to-date proposal-environment
    prepare_pbe()

    ### making variables and parameters available within the function
    le <- environment()
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v

    ### change to proposal state
    p$wh.slope <- exp(log(p$wh.slope) + rnorm(1, 0, h$si.wh))

    ### update proposal environment
    copy2pbe1("p", le)

    ### calculate proposalratio
    logproposalratio <- log(p$wh.slope) - log(pbe0$p$wh.slope)

    ### calculate likelihood
    propose_pbe("wh.slope")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal + logproposalratio + 
      dgamma(pbe1$p$wh.slope, shape = h$wh.s.sh, scale = h$wh.s.av/h$wh.s.sh, log = TRUE) - 
      dgamma(pbe0$p$wh.slope, shape = h$wh.s.sh, scale = h$wh.s.av/h$wh.s.sh, log = TRUE)
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
        accept_pbe("wh.slope")
    }
}

update_wh_exponent <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  ### change to proposal state
  p$wh.exponent <- exp(log(p$wh.exponent) + rnorm(1, 0, h$si.wh))
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$wh.exponent) - log(pbe0$p$wh.exponent)
  
  ### calculate likelihood
  propose_pbe("wh.exponent")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal + logproposalratio + 
    dgamma(pbe1$p$wh.exponent, shape = h$wh.e.sh, scale = h$wh.e.av/h$wh.e.sh, log = TRUE) - 
    dgamma(pbe0$p$wh.exponent, shape = h$wh.e.sh, scale = h$wh.e.av/h$wh.e.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("wh.exponent")
  }
}

update_wh_level <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  ### change to proposal state
  p$wh.level <- exp(log(p$wh.level) + rnorm(1, 0, h$si.wh))
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$wh.level) - log(pbe0$p$wh.level)
  
  ### calculate likelihood
  propose_pbe("wh.level")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikcoal - pbe0$logLikcoal + logproposalratio + 
    dgamma(pbe1$p$wh.level, shape = h$wh.0.sh, scale = h$wh.0.av/h$wh.0.sh, log = TRUE) - 
    dgamma(pbe0$p$wh.level, shape = h$wh.0.sh, scale = h$wh.0.av/h$wh.0.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("wh.level")
  }
}

update_dist_exponent <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  ### change to proposal state
  p$dist.exponent <- 1 + exp(log(p$dist.exponent - 1) + rnorm(1, 0, h$si.dist))
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$dist.exponent - 1) - log(pbe0$p$dist.exponent - 1)
  
  ### calculate likelihood
  propose_pbe("dist.exponent")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
    dgamma(pbe1$p$dist.exponent - 1, shape = h$dist.e.sh, scale = h$dist.e.av/h$dist.e.sh, log = TRUE) - 
    dgamma(pbe0$p$dist.exponent - 1, shape = h$dist.e.sh, scale = h$dist.e.av/h$dist.e.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("dist.exponent")
  }
}

update_dist_scale <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  ### change to proposal state
  p$dist.scale <- exp(log(p$dist.scale) + rnorm(1, 0, h$si.dist))
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$dist.scale) - log(pbe0$p$dist.scale)
  
  ### calculate likelihood
  propose_pbe("dist.scale")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
    dgamma(pbe1$p$dist.scale, shape = h$dist.s.sh, scale = h$dist.s.av/h$dist.s.sh, log = TRUE) - 
    dgamma(pbe0$p$dist.scale, shape = h$dist.s.sh, scale = h$dist.s.av/h$dist.s.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("dist.scale")
  }
}

update_dist_mean <- function() {
  ### create an up-to-date proposal-environment
  prepare_pbe()
  
  ### making variables and parameters available within the function
  le <- environment()
  h <- pbe0$h
  p <- pbe1$p
  v <- pbe1$v
  
  ### change to proposal state
  p$dist.mean <- exp(log(p$dist.mean) + rnorm(1, 0, h$si.dist))
  
  ### update proposal environment
  copy2pbe1("p", le)
  
  ### calculate proposalratio
  logproposalratio <- log(p$dist.mean) - log(pbe0$p$dist.mean)
  
  ### calculate likelihood
  propose_pbe("dist.mean")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikdist - pbe0$logLikdist + logproposalratio + 
    dgamma(pbe1$p$dist.mean, shape = h$dist.m.sh, scale = h$dist.m.av/h$dist.m.sh, log = TRUE) - 
    dgamma(pbe0$p$dist.mean, shape = h$dist.m.sh, scale = h$dist.m.av/h$dist.m.sh, log = TRUE)
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("dist.mean")
  }
}
