#The environments .phybreakenv and .phybreakenv.prop will contain sets
#of state variables + likelihood calculations. .phybreakenv is consistent
#with the current state of updating phybreak.object; .phybreakenv.prop
#is used to propose state and calculate likelihoods for these proposals.

#likarray is an array with likelihoods per nucleotide per SNP
#per node (sampling and coalescence). They are stored as vectors of
#length Nnodes*nSNPs*4, with implicit dim = c(4,nSNPs,Nnodes).

#Further, the environments contain SNP and SNPfreqs which
#are the underlying SNP data, and obs the outbreak size. Variable states
#are in the vectors nodetimes and nodeparents, in mutation rate mu,
#and in the calculated logLik. These are all consistent with likarray.
#logLikgen and logLiksam are the log-likelihood values of generation
#and sampling times with means mG and mS. The
#function '.build.phybreakenv' will be used to initialise, the function
# .propose.phybreakenv' will be used to change .phybreakenv.prop based
#on .phybreakenv and a proposed phybreak.object. The function
# .accept.phybreakenv' is to change .phybreakenv by copying .phybreakenv.prop.
#Finally, the functions 'changemS.phybreakenv' and 'changemG.phybreakenv'
#change the values of mS/mG and logLikmS/logLikmG when these have been resampled.

#fixed parameters
tinf.prop.shape <- 2
whrate0.prior.exp <- 1
whslope.prior.exp <- 1

#The environments are only used during MCMC-updating.

.phybreakenv <- new.env()
.phybreakenv.prop <- new.env()


.build.phybreakenv <- function(phybreak.obj) {
  .phybreakenv.prop$d <- phybreak.obj$d
  .phybreakenv.prop$v <- phybreak.obj$v
  .phybreakenv.prop$p <- phybreak.obj$p
  .phybreakenv.prop$h <- phybreak.obj$h
  .phybreakenv$tinf.prop.shape <- tinf.prop.shape
  .phybreakenv$d <- phybreak.obj$d
  .phybreakenv$v <- phybreak.obj$v
  .phybreakenv$p <- phybreak.obj$p
  .phybreakenv$h <- phybreak.obj$h


  .phybreakenv.prop$likarray <- array(1, dim = c(4,length(.phybreakenv.prop$d$SNPfr),
                                           2*.phybreakenv.prop$p$obs-1))
  .phybreakenv.prop$likarray[cbind(1,
                             rep(1:length(.phybreakenv.prop$d$SNPfr),
                                 .phybreakenv.prop$p$obs),
                             rep(1:.phybreakenv.prop$p$obs,
                                 each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "a" | .phybreakenv.prop$d$SNP == "n")
  .phybreakenv.prop$likarray[cbind(2,
                             rep(1:length(.phybreakenv.prop$d$SNPfr),
                                 .phybreakenv.prop$p$obs),
                             rep(1:.phybreakenv.prop$p$obs,
                                 each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "c" | .phybreakenv.prop$d$SNP == "n")
  .phybreakenv.prop$likarray[cbind(3,
                             rep(1:length(.phybreakenv.prop$d$SNPfr),
                                 .phybreakenv.prop$p$obs),
                             rep(1:.phybreakenv.prop$p$obs,
                                 each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "g" | .phybreakenv.prop$d$SNP == "n")
  .phybreakenv.prop$likarray[cbind(4,
                                   rep(1:length(.phybreakenv.prop$d$SNPfr),
                                       .phybreakenv.prop$p$obs),
                                   rep(1:.phybreakenv.prop$p$obs,
                                       each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "t" | .phybreakenv.prop$d$SNP == "n")


  .likseqenv(.phybreakenv.prop,
             (.phybreakenv.prop$p$obs + 1):(2*.phybreakenv.prop$p$obs - 1),
            1:.phybreakenv.prop$p$obs) #calculate logLik and adjust likarray

  .phybreakenv.prop$logLiksam <- with(phybreak.obj,
                                     .lik.sampletimes(p$shape.sample,
                                                      p$mean.sample,
                                                      v$nodetimes,
                                                      v$nodetypes
                                     ))

  .phybreakenv.prop$logLikgen <- with(phybreak.obj,
                                      .lik.gentimes(p$obs, p$shape.gen,
                                                    p$mean.gen, v$nodetimes,
                                                    v$nodehosts, v$nodetypes
                                      ))

  .phybreakenv.prop$logLikcoal <- with(phybreak.obj,
                                       .lik.coaltimes(p$obs, p$wh.model, p$wh.rate0, p$wh.lambda,
                                                      p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes
                                      ))

  .accept.phybreakenv()
}

.destroy.phybreakenv <- function(phybreak.obj.samples) {
  res <- list(
    d = .phybreakenv$d,
    v = .phybreakenv$v,
    p = .phybreakenv$p,
    h = .phybreakenv$h,
    s = phybreak.obj.samples
  )
  rm(list=ls(.phybreakenv), envir = .phybreakenv)
  rm(list=ls(.phybreakenv.prop), envir = .phybreakenv.prop)
  return(res)
}


.accept.phybreakenv <- function() {
  .phybreakenv$d <- .phybreakenv.prop$d
  .phybreakenv$v <- .phybreakenv.prop$v
  .phybreakenv$p <- .phybreakenv.prop$p
  .phybreakenv$h <- .phybreakenv.prop$h
  .phybreakenv$likarray <- .phybreakenv.prop$likarray
  .phybreakenv$logLik <- .phybreakenv.prop$logLik
  .phybreakenv$logLiksam <- .phybreakenv.prop$logLiksam
  .phybreakenv$logLikgen <- .phybreakenv.prop$logLikgen
  .phybreakenv$logLikcoal <- .phybreakenv.prop$logLikcoal
}

.prepare.phybreakenv <- function() {
  .phybreakenv.prop$d <- .phybreakenv$d
  .phybreakenv.prop$v <- .phybreakenv$v
  .phybreakenv.prop$p <- .phybreakenv$p
  .phybreakenv.prop$h <- .phybreakenv$h
  .phybreakenv.prop$likarray <- .phybreakenv$likarray
  .phybreakenv.prop$logLik <- .phybreakenv$logLik
  .phybreakenv.prop$logLiksam <- .phybreakenv$logLiksam
  .phybreakenv.prop$logLikgen <- .phybreakenv$logLikgen
  .phybreakenv.prop$logLikcoal <- .phybreakenv$logLikcoal
}

.propose.phybreakenv <- function() {
  if(.phybreakenv.prop$p$mu == .phybreakenv$p$mu) {
    #identify changed nodes
    chnodes <- which((.phybreakenv.prop$v$nodeparents != .phybreakenv$v$nodeparents) |
                       (.phybreakenv.prop$v$nodetimes != .phybreakenv$v$nodetimes))
    chnodes <- unique(unlist(sapply(chnodes,.ptr,
                                    pars=.phybreakenv.prop$v$nodeparents)))
    chnodes <- chnodes[chnodes > .phybreakenv.prop$p$obs &
                         chnodes < 2*.phybreakenv.prop$p$obs]
    #identify nodetips
    nodetips <- c(match(chnodes,.phybreakenv.prop$v$nodeparents),
                  3*.phybreakenv.prop$p$obs-
                    match(chnodes,rev(.phybreakenv.prop$v$nodeparents)))
    nodetips[nodetips >= 2*.phybreakenv.prop$p$obs] <-
      match(nodetips[nodetips >= 2*.phybreakenv.prop$p$obs],
            .phybreakenv.prop$v$nodeparents)
    nodetips <- nodetips[is.na(match(nodetips,chnodes))]
  } else {
    chnodes <- (.phybreakenv.prop$p$obs + 1):(2*.phybreakenv.prop$p$obs - 1)
    nodetips <- 1:.phybreakenv.prop$p$obs
  }

  .phybreakenv.prop$likarray <- .phybreakenv$likarray + 0

  if(!is.null(chnodes)) {
    .likseqenv(.phybreakenv.prop, chnodes, nodetips)
  } else {
    .phybreakenv.prop$logLik <- .phybreakenv$logLik
  }

  if(!all(.phybreakenv.prop$v$nodetimes == .phybreakenv$v$nodetimes) |
     !all(.phybreakenv.prop$v$nodeparents == .phybreakenv$v$nodeparents) |
     .phybreakenv.prop$p$mean.gen != .phybreakenv$p$mean.gen) {
    .phybreakenv.prop$logLikgen <- .lik.gentimes(
      .phybreakenv.prop$p$obs, .phybreakenv.prop$p$shape.gen,
      .phybreakenv.prop$p$mean.gen, .phybreakenv.prop$v$nodetimes,
      .phybreakenv.prop$v$nodehosts, .phybreakenv.prop$v$nodetypes
    )
  } else {
    .phybreakenv.prop$logLikgen <- .phybreakenv$logLikgen
  }

  if(!all(.phybreakenv.prop$v$nodetimes == .phybreakenv$v$nodetimes) |
     !all(.phybreakenv.prop$v$nodeparents == .phybreakenv$v$nodeparents) |
     .phybreakenv.prop$p$mean.sample != .phybreakenv$p$mean.sample) {
    .phybreakenv.prop$logLiksam <- .lik.sampletimes(
      .phybreakenv.prop$p$shape.sample, .phybreakenv.prop$p$mean.sample,
      .phybreakenv.prop$v$nodetimes, .phybreakenv.prop$v$nodetypes
    )
  } else {
    .phybreakenv.prop$logLiksam <- .phybreakenv$logLiksam
  }
  
  if(!all(.phybreakenv.prop$v$nodetimes == .phybreakenv$v$nodetimes) |
     !all(.phybreakenv.prop$v$nodeparents == .phybreakenv$v$nodeparents) |
     .phybreakenv.prop$p$wh.rate0 != .phybreakenv$p$wh.rate0 |
     .phybreakenv.prop$p$wh.slope != .phybreakenv$p$wh.slope) {
    .phybreakenv.prop$logLikcoal <- .lik.coaltimes(
      .phybreakenv.prop$p$obs, .phybreakenv.prop$p$wh.model, 
      .phybreakenv.prop$p$wh.rate0, .phybreakenv.prop$p$wh.lambda,
      .phybreakenv.prop$p$wh.slope, .phybreakenv.prop$v$nodetimes, 
      .phybreakenv.prop$v$nodehosts, .phybreakenv.prop$v$nodetypes
    )
  } else {
    .phybreakenv.prop$logLikcoal <- .phybreakenv$logLikcoal
  }
  

}



logLik.phybreak <- function(phybreak.obj, seq = TRUE, t.gen = FALSE,
                            t.sample = FALSE, t.coal = FALSE) {
  res <- 0
  if(seq) {
    res <- res + with(phybreak.obj, .likseq(t(d$SNP), d$SNPfr,
         v$nodeparents, v$nodetimes, p$mu,p$obs))
  } 
  if(t.gen) {
    res <- res + with(phybreak.obj, .lik.gentimes(p$obs, p$shape.gen, p$mean.gen,
         v$nodetimes, v$nodehosts, v$nodetypes))
  } 
  if(t.sample) {
    res <- res + with(phybreak.obj, .lik.sampletimes(p$shape.sample, p$mean.sample,
         v$nodetimes, v$nodetypes))
  } 
  if(t.coal) {
    res <- res + with(phybreak.obj, .lik.coaltimes(p$obs, p$wh.model, p$wh.rate0, p$wh.lambda,
                                      p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes))
  }
  return(res)
}


###For updating the global outbreak variables
burnin.phybreak <- function(phybreak.object, nburnin) {
#   parked.samples <- phybreak.object$s
#   res <- phybreak.object
#   res$s <- c()
  if(is.null(phybreak.object$h$est.wh)) {
    phybreak.object$h$si.wh <- if(wh.model == 3) {
      c(rep(c(wh.rate0,2*wh.rate0),50),rep(NA,900))
      } else {
        c(rep(c(wh.slope,2*wh.slope),50),rep(NA,900))
      }
    if(is.null(phybreak.object$h$est.r0)) {
      phybreak.object$h$est.wh <- FALSE
      phybreak.object$s$r0 <- c()
      phybreak.object$s$slope <- c()
    } else {
      phybreak.object$h$est.wh <- phybreak.object$h$est.r0
      phybreak.object$s$slope <- c()
    }
  }
  
  .build.phybreakenv(phybreak.object)


  for(rep in 1:nburnin) {
    for(i in sample(phybreak.object$p$obs)) {
      .updatehost(i)
    }
    if(phybreak.object$h$est.mG) .update.mG()
    if(phybreak.object$h$est.mS) .update.mS()
    if(phybreak.object$h$est.wh) .update.wh()
    .update.mu()
  }

  res <- .destroy.phybreakenv(phybreak.object$s)



  return(res)
}

sample.phybreak <- function(phybreak.object, nsample, thin) {
  if(is.null(phybreak.object$h$est.wh)) {
    phybreak.object$h$si.wh <- if(wh.model == 3) {
      c(rep(c(wh.rate0,2*wh.rate0),50),rep(NA,900))
    } else {
      c(rep(c(wh.slope,2*wh.slope),50),rep(NA,900))
    }
    if(is.null(phybreak.object$h$est.r0)) {
      phybreak.object$h$est.wh <- FALSE
      phybreak.object$s$r0 <- c()
      phybreak.object$s$slope <- c()
    } else {
      phybreak.object$h$est.wh <- phybreak.object$h$est.r0
      phybreak.object$s$slope <- c()
    }
  }
  
  s.post <- list(
    nodetimes = with(phybreak.object,cbind(s$nodetimes, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodehosts = with(phybreak.object,cbind(s$nodehosts, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodeparents = with(phybreak.object,cbind(s$nodeparents, matrix(NA, nrow=3*p$obs - 1, ncol=nsample))),
    mu = c(phybreak.object$s$mu, rep(NA, nsample)),
    mG = c(phybreak.object$s$mG, rep(NA, nsample)),
    mS = c(phybreak.object$s$mS, rep(NA, nsample)),
    r0 = c(phybreak.object$s$r0, rep(NA, nsample)),
    slope = c(phybreak.object$s$slope, rep(NA, nsample)),
    logLik = c(phybreak.object$s$logLik, rep(NA, nsample))
  )

  .build.phybreakenv(phybreak.object)

  for(sa in tail(1:length(s.post$mu), nsample)) {
    for(rep in 1:thin) {
      for(i in sample(phybreak.object$p$obs)) {
        .updatehost(i)
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
    s.post$r0[sa] <- .phybreakenv$p$wh.rate0
    s.post$slope[sa] <- .phybreakenv$p$wh.slope
    s.post$logLik[sa] <- .phybreakenv$logLik + .phybreakenv$logLiksam + 
      + .phybreakenv$logLikgen + .phybreakenv$logLikcoal
  }

  res <- .destroy.phybreakenv(s.post)


  return(res)

}


burnin.phybreak.plus <- function(phybreak.object, nburnin) {
  #   parked.samples <- phybreak.object$s
  #   res <- phybreak.object
  #   res$s <- c()
  if(is.null(phybreak.object$h$est.wh)) {
    phybreak.object$h$si.wh <- if(wh.model == 3) {
      c(rep(c(wh.rate0,2*wh.rate0),50),rep(NA,900))
    } else {
      c(rep(c(wh.slope,2*wh.slope),50),rep(NA,900))
    }
    if(is.null(phybreak.object$h$est.r0)) {
      phybreak.object$h$est.wh <- FALSE
      phybreak.object$s$r0 <- c()
      phybreak.object$s$slope <- c()
    } else {
      phybreak.object$h$est.wh <- phybreak.object$h$est.r0
      phybreak.object$s$slope <- c()
    }
  }
  
  .build.phybreakenv(phybreak.object)
  
  
  for(rep in 1:nburnin) {
    for(i in sample(phybreak.object$p$obs)) {
      if(runif(1) < 0.8) .updatehost(i) else .updatehost.keepphylo(i)
    }
    if(phybreak.object$h$est.mG) .update.mG()
    if(phybreak.object$h$est.mS) .update.mS()
    if(phybreak.object$h$est.wh) .update.wh()
    .update.mu()
  }
  
  res <- .destroy.phybreakenv(phybreak.object$s)
  
  
  
  return(res)
}

sample.phybreak.plus <- function(phybreak.object, nsample, thin) {
  if(is.null(phybreak.object$h$est.wh)) {
    phybreak.object$h$si.wh <- if(wh.model == 3) {
      c(rep(c(wh.rate0,2*wh.rate0),50),rep(NA,900))
    } else {
      c(rep(c(wh.slope,2*wh.slope),50),rep(NA,900))
    }
    if(is.null(phybreak.object$h$est.r0)) {
      phybreak.object$h$est.wh <- FALSE
      phybreak.object$s$r0 <- c()
      phybreak.object$s$slope <- c()
    } else {
      phybreak.object$h$est.wh <- phybreak.object$h$est.r0
      phybreak.object$s$slope <- c()
    }
  }
  
  s.post <- list(
    nodetimes = with(phybreak.object,cbind(s$nodetimes, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodehosts = with(phybreak.object,cbind(s$nodehosts, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodeparents = with(phybreak.object,cbind(s$nodeparents, matrix(NA, nrow=3*p$obs - 1, ncol=nsample))),
    mu = c(phybreak.object$s$mu, rep(NA, nsample)),
    mG = c(phybreak.object$s$mG, rep(NA, nsample)),
    mS = c(phybreak.object$s$mS, rep(NA, nsample)),
    r0 = c(phybreak.object$s$r0, rep(NA, nsample)),
    slope = c(phybreak.object$s$slope, rep(NA, nsample)),
    logLik = c(phybreak.object$s$logLik, rep(NA, nsample))
  )
  
  .build.phybreakenv(phybreak.object)
  
  for(sa in tail(1:length(s.post$mu), nsample)) {
    for(rep in 1:thin) {
      for(i in sample(phybreak.object$p$obs)) {
        if(runif(1) < .8) .updatehost(i) else .updatehost.keepphylo(i)
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
    s.post$r0[sa] <- .phybreakenv$p$wh.rate0
    s.post$slope[sa] <- .phybreakenv$p$wh.slope
    s.post$logLik[sa] <- .phybreakenv$logLik + .phybreakenv$logLiksam + 
      + .phybreakenv$logLikgen + .phybreakenv$logLikcoal
  }
  
  res <- .destroy.phybreakenv(s.post)
  
  
  return(res)
  
}

burnin.phybreak.alt <- function(phybreak.object, nburnin) {
  #   parked.samples <- phybreak.object$s
  #   res <- phybreak.object
  #   res$s <- c()
  if(is.null(phybreak.object$h$est.wh)) {
    phybreak.object$h$si.wh <- if(wh.model == 3) {
      c(rep(c(wh.rate0,2*wh.rate0),50),rep(NA,900))
    } else {
      c(rep(c(wh.slope,2*wh.slope),50),rep(NA,900))
    }
    if(is.null(phybreak.object$h$est.r0)) {
      phybreak.object$h$est.wh <- FALSE
      phybreak.object$s$r0 <- c()
      phybreak.object$s$slope <- c()
    } else {
      phybreak.object$h$est.wh <- phybreak.object$h$est.r0
      phybreak.object$s$slope <- c()
    }
  }
  
  .build.phybreakenv(phybreak.object)
  
  
  for(rep in 1:nburnin) {
    for(i in sample(phybreak.object$p$obs)) {
      if(runif(1) < 0.8) .updatehost(i) else .updatehost.keepphylo.alt(i)
    }
    if(phybreak.object$h$est.mG) .update.mG()
    if(phybreak.object$h$est.mS) .update.mS()
    if(phybreak.object$h$est.wh) .update.wh()
    .update.mu()
  }
  
  res <- .destroy.phybreakenv(phybreak.object$s)
  
  
  
  return(res)
}

sample.phybreak.alt <- function(phybreak.object, nsample, thin) {
  if(is.null(phybreak.object$h$est.wh)) {
    phybreak.object$h$si.wh <- if(wh.model == 3) {
      c(rep(c(wh.rate0,2*wh.rate0),50),rep(NA,900))
    } else {
      c(rep(c(wh.slope,2*wh.slope),50),rep(NA,900))
    }
    if(is.null(phybreak.object$h$est.r0)) {
      phybreak.object$h$est.wh <- FALSE
      phybreak.object$s$r0 <- c()
      phybreak.object$s$slope <- c()
    } else {
      phybreak.object$h$est.wh <- phybreak.object$h$est.r0
      phybreak.object$s$slope <- c()
    }
  }
  
  s.post <- list(
    nodetimes = with(phybreak.object,cbind(s$nodetimes, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodehosts = with(phybreak.object,cbind(s$nodehosts, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodeparents = with(phybreak.object,cbind(s$nodeparents, matrix(NA, nrow=3*p$obs - 1, ncol=nsample))),
    mu = c(phybreak.object$s$mu, rep(NA, nsample)),
    mG = c(phybreak.object$s$mG, rep(NA, nsample)),
    mS = c(phybreak.object$s$mS, rep(NA, nsample)),
    r0 = c(phybreak.object$s$r0, rep(NA, nsample)),
    slope = c(phybreak.object$s$slope, rep(NA, nsample)),
    logLik = c(phybreak.object$s$logLik, rep(NA, nsample))
  )
  
  .build.phybreakenv(phybreak.object)
  
  for(sa in tail(1:length(s.post$mu), nsample)) {
    for(rep in 1:thin) {
      for(i in sample(phybreak.object$p$obs)) {
        if(runif(1) < .8) .updatehost(i) else .updatehost.keepphylo.alt(i)
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
    s.post$r0[sa] <- .phybreakenv$p$wh.rate0
    s.post$slope[sa] <- .phybreakenv$p$wh.slope
    s.post$logLik[sa] <- .phybreakenv$logLik + .phybreakenv$logLiksam + 
      + .phybreakenv$logLikgen + .phybreakenv$logLikcoal
  }
  
  res <- .destroy.phybreakenv(s.post)
  
  
  return(res)
  
}



testupdate <- function(phybreak.object, hostID) {
  .build.phybreakenv(phybreak.object)
  .updatehost.keepphylo(hostID)
  res <- .destroy.phybreakenv(phybreak.object$s)
  return(res)
  
}


.updatehost <- function(hostID) {
  .prepare.phybreakenv()
  .phybreakenv.prop$hostID <- hostID
  evalq({
    #    tinf.prop <- v$nodetimes[hostID] - rgamma(1,shape=p$shape.sample,scale=p$mean.sample/p$shape.sample)
    tinf.prop <- v$nodetimes[hostID] - 
      rgamma(1,shape=tinf.prop.shape,scale=p$mean.sample/tinf.prop.shape)

    ##going down the decision tree
    if(tinf.prop < suppressWarnings(min(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID]))) {
      #Q1=Y : tinf.prop before first infectee
      if(v$nodehosts[hostID + 2*p$obs - 1] == 0) {
        #Q12=YY : tinf.prop before first infectee & hostID is index
#        print(paste0(hostID,": path A"))
        .updatepathA()
      } else {
        #Q12=YN : tinf.prop before first infectee & hostID is not index
        if(tinf.prop < v$nodetimes[v$nodehosts == 0]) {
          #Q123=YNY : hostID is not index & tinf.prop is before infection of index
#          print(paste0(hostID,": path B"))
          .updatepathB()
#           iees <- which(v$nodehosts == (which(v$nodehosts == 0) - 2*p$obs + 1) & v$nodetypes == "t")
#           if(hostID == iees[order(v$nodetimes[iees])][1] - 2*p$obs + 1) {
#             #Q1234=YNYY : tinf.prop before infection of index & hostID is index's first infectee
#             .updatepathB()
#           } #Q1234=YNYN : tinf.prop before infection of index & hostID is not index nor index's first infectee
#             ### DO NOTHING
        } else {
          #Q123=YNN : tinf.prop before first infectee and after infection of index
#          print(paste0(hostID,": path C"))
          .updatepathC()
        }
      }
    } else
    {
      #Q1=N : tinf.prop after first infectee
      if(v$nodehosts[hostID + 2*p$obs - 1] == 0) {
        #Q12=NY : tinf.prop after first infectee & hostID is index
        if(tinf.prop <= sort(c(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID],Inf))[2]) {
#           if(tinf.prop <= sort(tail(v$nodetimes,p$obs))[3]) {
#            #Q123=NYY : hostID is index & tinf.prop is after first infectee but before third infection in outbreak
            #Q123=NYY : hostID is index & tinf.prop is after first infectee but before second infectee
#          print(paste0(hostID,": path D"))
          .updatepathD()
        } else {
          #Q123=NYN : hostID is index & tinf.prop is after third infection in outbreak
#          print(paste0(hostID,": path E1"))
          if(runif(1,0,1) < .5) .updatepathE() else .updatepathF()
        }
      } else {
        #Q12=NN : tinf.prop after first infectee & hostID is not index
#        print(paste0(hostID,": path E2"))
        if(runif(1,0,1) < .5) .updatepathE() else .updatepathF()
      }
    }
  }, .phybreakenv.prop)
}

.updatehost.keepphylo <- function(hostID) {
  .prepare.phybreakenv()
  .phybreakenv.prop$hostID <- hostID
  evalq({
    mean.propose <- p$mean.sample
    tinf.prop <- v$nodetimes[hostID] - 
      rgamma(1,shape=tinf.prop.shape,scale=mean.propose/tinf.prop.shape)
    hostiorID <- v$nodehosts[2*p$obs - 1 + hostID]
    
#     print(paste0("hostID: ",hostID))
#     print(paste0("hostiorID: ",hostiorID))
#     print(paste0("tinf.prop = ",tinf.prop))
    
    if(hostiorID == 0) {
      #hostID is index case, test if tinf.prop is before first coalescent node
#       print(paste0("time first coal node= ",v$nodetimes[which(v$nodeparents == 2*p$obs - 1 + hostID)]))
      if(tinf.prop < v$nodetimes[which(v$nodeparents == 2*p$obs - 1 + hostID)]) {
#         print("take path A")
        .updatepathA.phylo()
      } #else reject
    } else {
      hostioriorID <- v$nodehosts[2*p$obs - 1 + hostiorID]
      nodemrca <- .ptr(v$nodeparents, hostID)[.ptr(v$nodeparents, hostID)%in%.ptr(v$nodeparents,hostiorID)][1]
      timemrca <- v$nodetimes[nodemrca]
#       print(paste0("hostioriorID: ",hostioriorID))
#       print(paste0("timemrca = ",timemrca))
      if(hostioriorID == 0) {
        #hostID's infector is index case, test if tinf.prop is before first coalescent node
#         print(paste0("time first coal node= ",v$nodetimes[which(v$nodeparents == 2*p$obs - 1 + hostiorID)]))
                if(tinf.prop < timemrca) {
          tinf2.prop <- v$nodetimes[hostiorID] - 
            rgamma(1,shape=tinf.prop.shape,scale=p$mean.sample/tinf.prop.shape)
#           print(paste0("tinf2.prop = ",tinf2.prop))
          #test if tinf2.prop is after timemrca -> to have a proper transmission tree
          if(tinf2.prop > timemrca) {
            ##TO PROGRAMME: first change tinf.hostiorID, then change tinf.hostID, then test within-host likelihood##
#             print("take path D")
            tinf.prop <- v$nodetimes[hostiorID + 2*p$obs - 1]
            .updatepathD.phylo()
          } #else reject
        } else if(tinf.prop > timemrca) {
          ##TO PROGRAMME: change tinf.hostID, then test within-host likelihood##
#           print("take path C")
          .updatepathC.phylo()
        } #else reject
      } else {
        #hostID's infector is not index case, determine mrca with infector's infector
        nodemrcaior <- .ptr(v$nodeparents, hostID)[.ptr(v$nodeparents, hostID)%in%.ptr(v$nodeparents,hostioriorID)][1]
        timemrcaior <- v$nodetimes[nodemrcaior]
#         print(paste0("timemrcaior = ",timemrcaior))
        
        #test if tinf.prop is after timemrcaior
        if(tinf.prop > timemrcaior) {
          #test if tinf.prop is before timemrca
          if(tinf.prop < timemrca) {
            tinf2.prop <- v$nodetimes[hostiorID] - 
              rgamma(1,shape=tinf.prop.shape,scale=p$mean.sample/tinf.prop.shape)
#             print(paste0("tinf2.prop = ",tinf2.prop))
            #test if tinf2.prop is after timemrca -> to have a proper transmission tree
            if(tinf2.prop > timemrca) {
              ##TO PROGRAMME: first change tinf.hostiorID, then change tinf.hostID, then test within-host likelihood##
#               print("take path B")
              .updatepathB.phylo()
            } #else reject
          } else {
            ##TO PROGRAMME: change tinf.hostID, then test within-host likelihood##
#             print("take path C")
            .updatepathC.phylo()
          }
        } #else reject
      }
    }
    

  }, .phybreakenv.prop)
  
  
}

.updatehost.keepphylo.alt <- function(hostID) {
  .prepare.phybreakenv()
  .phybreakenv.prop$hostID <- hostID
  evalq({
    nodeUC <- v$nodeparents[2*p$obs - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2*p$obs - 1 + hostID)
    hostiorID <- v$nodehosts[2*p$obs - 1 + hostID]
    if(nodeUC != 0) {
      intervalNOW <- v$nodetimes[nodeDC] - v$nodetimes[nodeUC]
      nodeUUC <- v$nodeparents[nodeUC]
      intervalUP <- v$nodetimes[nodeUC] - v$nodetimes[nodeUUC]
    } else nodeUUC <- -1
    if(any(.ptr(v$nodeparents,hostiorID)==nodeUC) & nodeUC != 0) {
      nodeUDC <- which(v$nodeparents == nodeUC)
      nodeUDC <- nodeUDC[nodeUDC != hostID + 2*p$obs - 1]
      intervalUPDOWN <- v$nodetimes[nodeUDC] - v$nodetimes[nodeUC]
      nodeUP <- v$nodeparents[hostiorID + 2*p$obs - 1]
      nodeDP <- which(v$nodeparents == 2*p$obs - 1 + hostiorID)
    } else nodeUDC <- -1
    if(nodeDC != hostID) {
      nodeDDC <- intersect(which(v$nodeparents == nodeDC),.ptr(v$nodeparents,hostID))
      intervalDOWN <- v$nodetimes[nodeDDC] - v$nodetimes[nodeDC]
    } else nodeDDC <- -1
    if(nodeUC == 0) {
      nodeDDC2 <- which(v$nodeparents == nodeDC)
      nodeDDC2 <- nodeDDC2[nodeDDC2 != nodeDDC]
      intervalNOW <- intervalUPDOWN <- 1
    }
    
    direction <- sample(1:3,1)
    
    if(direction == 1) {
      if(nodeUUC != -1) {
        if(nodeUDC == -1) {
          tinf.prop <- runif(1, v$nodetimes[nodeUUC], v$nodetimes[nodeUC])
          .updatepathUP()
        } else {
          tinf.prop <- runif(1, v$nodetimes[nodeUC], v$nodetimes[nodeUDC])
          .updatepathUPDOWN()
        }
      }
    } else if(direction == 2) {
      if(nodeDDC != -1) {
        if(nodeUUC != -1) {
          tinf.prop <- runif(1, v$nodetimes[nodeDC], v$nodetimes[nodeDDC])
          .updatepathDOWN()
        } else {
          if(nodeDDC2 >= 2*p$obs) {
            hostiorID <- hostID
            hostID <- nodeDDC2 + 1 - 2*p$obs
            nodeUP <- nodeUC
            nodeDP <- nodeDC
            nodeUC <- nodeDP
            nodeDC <- which(v$nodeparents == nodeDDC2)
            nodeUDC <- which(v$nodeparents == nodeUC)
            nodeUDC <- nodeUDC[nodeUDC != hostID + 2*p$obs - 1]
            tinf.prop <- runif(1, v$nodetimes[nodeUC], v$nodetimes[nodeUDC])
            .updatepathUPDOWN()
          }
        }
      }
    } else if(nodeUUC != -1) {
      tinf.prop <- runif(1, v$nodetimes[nodeUC], v$nodetimes[nodeDC])
      .updatepathNOW()
    }
    
    
    
  }, .phybreakenv.prop)
  
  
}


#helper functions for updating variables
{
.updatepathA <- function() {


  evalq({
    ###change to proposal state

    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop

    v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
      v$nodetimes[hostID + 2*p$obs - 1] +                      #...to the infection time +
      .samplecoaltimes(v$nodetimes[v$nodehosts == hostID & v$nodetypes != "c"] - v$nodetimes[hostID + 2*p$obs - 1],
                       p$wh.model, p$wh.lambda, p$wh.rate0, p$wh.slope)  #...sampled coalescence times

    v$nodeparents[v$nodehosts == hostID] <-     #change the parent nodes of all nodes in hostID...
      .sampletopology(which(v$nodehosts == hostID), v$nodetimes[v$nodehosts == hostID], v$nodetypes[v$nodehosts == hostID],
                      hostID + 2*p$obs - 1, p$wh.model)
    #...to a correct topology, randomized where possible



    ###calculate acceptance probability
    logproposalratio <-
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log=TRUE)

  },
  .phybreakenv.prop)






  .propose.phybreakenv()

  evalq(
  logaccprob <- logLik + logLikgen +
    logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
    .phybreakenv$logLiksam + logproposalratio,
  .phybreakenv.prop)



  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
#    print("accepted")
    .accept.phybreakenv()
  }

}

.updatepathB <- function() {

  evalq({
    ###change to proposal state
    index.current.ID <- which(v$nodehosts == 0) - 2*p$obs + 1
    infector.current.ID <- v$nodehosts[hostID + 2*p$obs - 1]

    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop

    v$nodehosts[hostID + 2*p$obs - 1] <- 0
    v$nodeparents[hostID + 2*p$obs - 1] <- 0
    v$nodehosts[index.current.ID + 2*p$obs - 1] <- hostID
    v$nodehosts[v$nodehosts == infector.current.ID & v$nodetypes == "c"][1] <- hostID

    for(ID in c(hostID,infector.current.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0, p$wh.slope)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible

    }


    ###calculate acceptance probability
    dens.infectorcurrent <-
      (dgamma(.phybreakenv$v$nodetimes[hostID + 2*p$obs - 1] -
                tail(.phybreakenv$v$nodetimes,p$obs),
              shape = p$shape.gen,
              scale = p$mean.gen/p$shape.gen)) +
      (.phybreakenv$v$nodetimes[hostID + 2*p$obs - 1] -
         tail(.phybreakenv$v$nodetimes,p$obs) > 0) /
      h$dist[hostID,]
    dens.infectorcurrent[hostID] <- 0
    
    logproposalratio <-
      log(dens.infectorcurrent[infector.current.ID] /
            (sum(dens.infectorcurrent))) +
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log=TRUE)

  },.phybreakenv.prop)




  .propose.phybreakenv()

  evalq(
    logaccprob <- logLik + logLikgen +
      logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam + logproposalratio,
    .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
#    print("accepted")
    .accept.phybreakenv()
  }

}

.updatepathC <- function() {


  evalq({
    ###change to proposal state
    infector.current.ID <- v$nodehosts[hostID + 2*p$obs - 1]
    dens.infectorproposal <-
      (dgamma(tinf.prop - tail(v$nodetimes,p$obs),
              shape = p$shape.gen,
              scale = p$mean.gen/p$shape.gen)) +
      (tinf.prop - tail(v$nodetimes,p$obs) > 0) / h$dist[hostID,]
    dens.infectorproposal[hostID] <- 0
    infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)



    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    v$nodehosts[hostID + 2*p$obs - 1] <- infector.proposed.ID
    v$nodehosts[v$nodehosts == infector.current.ID & v$nodetypes == "c"][1] <- infector.proposed.ID


    for(ID in c(hostID,infector.current.ID,infector.proposed.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0, p$wh.slope)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible
    }


  ###calculate acceptance probability
  dens.infectorcurrent <-
    (dgamma(.phybreakenv$v$nodetimes[hostID + 2*p$obs - 1] -
              tail(.phybreakenv$v$nodetimes,p$obs),
            shape = p$shape.gen,
            scale = p$mean.gen/p$shape.gen)) +
    (.phybreakenv$v$nodetimes[hostID + 2*p$obs - 1] -
       tail(.phybreakenv$v$nodetimes,p$obs) > 0) /
    h$dist[hostID,]
  dens.infectorcurrent[hostID] <- 0
  logproposalratio <-
    log(dens.infectorcurrent[infector.current.ID] *
          sum(dens.infectorproposal) /
          (dens.infectorproposal[infector.proposed.ID] *
             sum(dens.infectorcurrent))) +
    dgamma(.phybreakenv$v$nodetimes[hostID] -
             .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
           shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape,log=TRUE) -
    dgamma(v$nodetimes[hostID] -
             v$nodetimes[hostID + 2*p$obs - 1],
           shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log=TRUE)

  },.phybreakenv.prop)




  .propose.phybreakenv()

  evalq(
  logaccprob <- logLik + logLikgen +
    logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
    .phybreakenv$logLiksam + logproposalratio,
  .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
#    print("accepted")
    .accept.phybreakenv()
  }

}

.updatepathD <- function() {


  evalq({
    ###change to proposal state
    dens.infectorproposal <-
      (dgamma(tinf.prop - tail(v$nodetimes,p$obs),
              shape = p$shape.gen,
              scale = p$mean.gen/p$shape.gen)) +
      (tinf.prop - tail(v$nodetimes,p$obs) > 0) / h$dist[hostID,]
    dens.infectorproposal[hostID] <- 0
    infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
    
    iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
    infectee.current.ID <-
      iees.nodeIDs[v$nodetimes[iees.nodeIDs] ==
                     min(v$nodetimes[iees.nodeIDs])] - 2*p$obs + 1
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop

    v$nodehosts[infectee.current.ID + 2*p$obs - 1] <- 0
    v$nodeparents[infectee.current.ID + 2*p$obs - 1] <- 0
    v$nodehosts[hostID + 2*p$obs - 1] <- infector.proposed.ID
    v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infector.proposed.ID

    for(ID in c(hostID,infector.proposed.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0, p$wh.slope)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible

    }

  ###calculate acceptance probability
  logproposalratio <-
    log(sum(dens.infectorproposal) /
          (dens.infectorproposal[infector.proposed.ID])) +
    dgamma(.phybreakenv$v$nodetimes[hostID] -
             .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
           shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape,log=TRUE) -
    dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
           shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log=TRUE)

  },.phybreakenv.prop)





  .propose.phybreakenv()

  evalq(
  logaccprob <- logLik + logLikgen +
    logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
    .phybreakenv$logLiksam + logproposalratio,
  .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
#    print("accepted")
    .accept.phybreakenv()
  }

}

.updatepathE <- function() {

  evalq({
    ###change to proposal state
    iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
    infectee.current.ID <-
      iees.nodeIDs[v$nodetimes[iees.nodeIDs] ==
                     min(v$nodetimes[iees.nodeIDs])] - 2*p$obs + 1
#    print(v$nodetimes[hostID]-v$nodetimes[iees.nodeIDs])
    infector.current.ID <- v$nodehosts[hostID + 2*p$obs - 1]
    inftime.current.hostID <- v$nodetimes[hostID + 2*p$obs - 1]

    v$nodetimes[hostID + 2*p$obs - 1] <- v$nodetimes[infectee.current.ID + 2*p$obs - 1]
    v$nodetimes[infectee.current.ID + 2*p$obs - 1] <- inftime.current.hostID

    v$nodehosts[hostID + 2*p$obs - 1] <- infectee.current.ID
    v$nodehosts[infectee.current.ID + 2*p$obs - 1] <- infector.current.ID
    v$nodeparents[infectee.current.ID + 2*p$obs - 1] <- v$nodeparents[hostID + 2*p$obs - 1]
    v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infectee.current.ID

    for(ID in c(hostID,infectee.current.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0, p$wh.slope)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible

    }

  ###calculate acceptance probability
  if(infector.current.ID != 0) {
    ##probability to sample after FIRST infectee
    logproposalratio <-
      pgamma(v$nodetimes[infectee.current.ID] - v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE) -
      pgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE)
  } else {
    ##probability to sample after SECOND infectee
#    print(v$nodetimes[hostID]-.phybreakenv$v$nodetimes[iees.nodeIDs])
    
    sampleinterval.current.ID <- v$nodetimes[hostID] - sort(.phybreakenv$v$nodetimes[iees.nodeIDs])[2]
    iees.newindex <- which(v$nodehosts == infectee.current.ID & v$nodetypes == "t")
    sampleinterval.newindex <- v$nodetimes[infectee.current.ID] - sort(c(v$nodetimes[iees.newindex],Inf))[2]
    
    
    logproposalratio <-
      pgamma(sampleinterval.newindex,
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE) -
      pgamma(sampleinterval.current.ID,
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE)
#    print(paste0(sampleinterval.current.ID, " ",sampleinterval.newindex, " ",logproposalratio))
    
  }

  },.phybreakenv.prop)




  .propose.phybreakenv()

  evalq({
    logaccprob <- logLik + logLikgen +
      logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam + logproposalratio
  },.phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
#    print("accepted")
    .accept.phybreakenv()
  }


}

.updatepathF <- function() {
  
  evalq({
    ###change to proposal state
    iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t") 
    infectee.current.ID <-
      iees.nodeIDs[v$nodetimes[iees.nodeIDs] ==
                     min(v$nodetimes[iees.nodeIDs])] - 2*p$obs + 1
    infecteenodes.primary <- iees.nodeIDs[iees.nodeIDs != (infectee.current.ID + 2*p$obs - 1)]
    infecteenodes.secondary <- which(v$nodehosts == infectee.current.ID & v$nodetypes == "t")
    #    print(v$nodetimes[hostID]-v$nodetimes[iees.nodeIDs])
    infector.current.ID <- v$nodehosts[hostID + 2*p$obs - 1]
    inftime.current.hostID <- v$nodetimes[hostID + 2*p$obs - 1]
    coalnodes.current.ID <- (v$nodehosts == hostID & v$nodetypes == "c")
    
    v$nodetimes[hostID + 2*p$obs - 1] <- v$nodetimes[infectee.current.ID + 2*p$obs - 1]
    v$nodetimes[infectee.current.ID + 2*p$obs - 1] <- inftime.current.hostID
    
    v$nodehosts[hostID + 2*p$obs - 1] <- infectee.current.ID
    v$nodehosts[infectee.current.ID + 2*p$obs - 1] <- infector.current.ID
    v$nodehosts[infecteenodes.primary] <- infectee.current.ID
    v$nodehosts[infecteenodes.secondary] <- hostID
    v$nodeparents[infectee.current.ID + 2*p$obs - 1] <- v$nodeparents[hostID + 2*p$obs - 1]
    v$nodehosts[v$nodehosts == infectee.current.ID & v$nodetypes == "c"] <- hostID
    v$nodehosts[coalnodes.current.ID] <- infectee.current.ID
    
    for(ID in c(hostID,infectee.current.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0, p$wh.slope)  #...sampled coalescence times
      
      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible
      
    }
    
    ###calculate acceptance probability
    if(infector.current.ID != 0) {
      ##probability to sample after FIRST infectee
      logproposalratio <-
        pgamma(v$nodetimes[infectee.current.ID] - v$nodetimes[hostID + 2*p$obs - 1],
               shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE) -
        pgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
               shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE)
    } else {
      ##probability to sample after SECOND infectee
      #    print(v$nodetimes[hostID]-.phybreakenv$v$nodetimes[iees.nodeIDs])
      
      sampleinterval.current.ID <- v$nodetimes[hostID] - sort(.phybreakenv$v$nodetimes[iees.nodeIDs])[2]
      sampleinterval.newindex <- v$nodetimes[infectee.current.ID] - sort(.phybreakenv$v$nodetimes[iees.nodeIDs])[2]
      
      
      logproposalratio <-
        pgamma(sampleinterval.newindex,
               shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE) -
        pgamma(sampleinterval.current.ID,
               shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log.p = TRUE)
      #    print(paste0(sampleinterval.current.ID, " ",sampleinterval.newindex, " ",logproposalratio))
      
    }
    
  },.phybreakenv.prop)
  
  
  
  
  .propose.phybreakenv()
  
  evalq({
    logaccprob <- logLik + logLikgen +
      logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam + logproposalratio
  },.phybreakenv.prop)
  

  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}


}

{
.updatepathA.phylo <- function() {

  
  evalq({
    ###change to proposal state
    
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    

    ###calculate acceptance probability
    logproposalratio <- 
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape, log=TRUE)
#       print(paste0("logpropratio = ",logproposalratio))
  },
  .phybreakenv.prop)
 
  .propose.phybreakenv() 
 
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
#     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  

  
  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}

.updatepathB.phylo <- function() {
  evalq({
    ###change to proposal state
    
    #first the infector
    v$nodetimes[hostiorID + 2*p$obs - 1] <- tinf2.prop
    nodesinvolved <- which(v$nodehosts == hostiorID | v$nodehosts == hostioriorID)
    nodeUC <- v$nodeparents[2*p$obs - 1 + hostiorID]
    nodeDC <- which(v$nodeparents == 2*p$obs - 1 + hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    nodeUP <- .ptr(v$nodeparents,hostiorID)[v$nodetimes[.ptr(v$nodeparents,hostiorID)] < tinf2.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP),.ptr(v$nodeparents,hostiorID))
    v$nodeparents[nodeDP] <- 2*p$obs - 1 + hostiorID
    v$nodeparents[2*p$obs - 1 + hostiorID] <- nodeUP
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostiorID)) {
        v$nodehosts[i] <- hostiorID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostiorID] <- hostioriorID
    
    #then the focal host itself
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostioriorID)
    nodeUC <- v$nodeparents[2*p$obs - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2*p$obs - 1 + hostID)
    v$nodeparents[nodeDC] <- nodeUC
    nodeUP <- c(.ptr(v$nodeparents,hostID)[v$nodetimes[.ptr(v$nodeparents,hostID)] < tinf.prop],0)[1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP),.ptr(v$nodeparents,hostID))
    v$nodeparents[nodeDP] <- 2*p$obs - 1 + hostID
    v$nodeparents[2*p$obs - 1 + hostID] <- nodeUP
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostID] <- hostioriorID
    
    
    
    ###calculate acceptance probability
    logproposalratio <- 
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape, log=TRUE) +
      dgamma(.phybreakenv$v$nodetimes[hostiorID] -
               .phybreakenv$v$nodetimes[hostiorID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape,log=TRUE) -
      dgamma(v$nodetimes[hostiorID] -
               v$nodetimes[hostiorID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log=TRUE)
#     print(paste0("logpropratio = ",logproposalratio))
    
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
#     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  

  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
}

.updatepathC.phylo <- function() {
  evalq({
    ###change to proposal state
    

    #only for the focal host itself
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostiorID)
    nodeUC <- v$nodeparents[2*p$obs - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2*p$obs - 1 + hostID)
    v$nodeparents[nodeDC] <- nodeUC
    nodeUP <- .ptr(v$nodeparents,hostID)[v$nodetimes[.ptr(v$nodeparents,hostID)] < tinf.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP),.ptr(v$nodeparents,hostID))
    v$nodeparents[nodeDP] <- 2*p$obs - 1 + hostID
    v$nodeparents[2*p$obs - 1 + hostID] <- nodeUP
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostiorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostID] <- hostiorID
    
    
    
    ###calculate acceptance probability
    logproposalratio <- 
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape, log=TRUE) 
#     print(paste0("logpropratio = ",logproposalratio))
    
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
#     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  

  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}

.updatepathD.phylo <- function() {
  evalq({

    ###change to proposal state
    
    #first the infector
    v$nodetimes[hostiorID + 2*p$obs - 1] <- tinf2.prop
    nodesinvolved <- which(v$nodehosts == hostiorID | v$nodehosts == hostioriorID)
    nodeUC <- v$nodeparents[2*p$obs - 1 + hostiorID]
    nodeDC <- which(v$nodeparents == 2*p$obs - 1 + hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    nodeUP <- .ptr(v$nodeparents,hostiorID)[v$nodetimes[.ptr(v$nodeparents,hostiorID)] < tinf2.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP),.ptr(v$nodeparents,hostiorID))
    v$nodeparents[nodeDP] <- 2*p$obs - 1 + hostiorID
    v$nodeparents[2*p$obs - 1 + hostiorID] <- nodeUP
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostiorID)) {
        v$nodehosts[i] <- hostiorID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostiorID] <- hostioriorID
    
    #then the focal host itself
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostioriorID)
    nodeUC <- v$nodeparents[2*p$obs - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2*p$obs - 1 + hostID)
    v$nodeparents[nodeDC] <- nodeUC
    nodeUP <- c(.ptr(v$nodeparents,hostID)[v$nodetimes[.ptr(v$nodeparents,hostID)] < tinf.prop],0)[1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP),.ptr(v$nodeparents,hostID))
    v$nodeparents[nodeDP] <- 2*p$obs - 1 + hostID
    v$nodeparents[2*p$obs - 1 + hostID] <- nodeUP
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostID] <- hostioriorID
    
    
    
    ###calculate acceptance probability
    logproposalratio <- 
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape,log=TRUE) -
      log(1 - pgamma(v$nodetimes[hostID] -
               timemrca,
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape)) +
      log(1 - pgamma(.phybreakenv$v$nodetimes[hostiorID] -
               timemrca,
             shape = tinf.prop.shape, scale = mean.propose/tinf.prop.shape)) -
      dgamma(v$nodetimes[hostiorID] -
               v$nodetimes[hostiorID + 2*p$obs - 1],
             shape = tinf.prop.shape, scale = p$mean.sample/tinf.prop.shape, log=TRUE)
#     print(paste0("logpropratio = ",logproposalratio))
    
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
#     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  
  
  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
}

}


{
.updatepathNOW <- function() {
  
  
  evalq({
    ###change to proposal state
    
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    
    
    ###calculate acceptance probability
    logproposalratio <- 0
    #       print(paste0("logpropratio = ",logproposalratio))
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
    #     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  
  
  
  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}


.updatepathUP <- function() {
  evalq({
    ###change to proposal state
    hostiorID <- v$nodehosts[hostID + 2*p$obs - 1]
    
    #only for the focal host itself
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    v$nodeparents[nodeUC] <- 2*p$obs - 1 + hostID
    v$nodeparents[2*p$obs - 1 + hostID] <- nodeUUC
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostiorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostID] <- hostiorID
    
    
    
    ###calculate acceptance probability
    logproposalratio <- log(intervalUP/intervalNOW)
    #     print(paste0("logpropratio = ",logproposalratio))
    
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
    #     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  
  
  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}

.updatepathDOWN <- function() {
  evalq({
    ###change to proposal state
    hostiorID <- v$nodehosts[hostID + 2*p$obs - 1]
    
    #only for the focal host itself
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    v$nodeparents[nodeDDC] <- 2*p$obs - 1 + hostID
    v$nodeparents[2*p$obs - 1 + hostID] <- nodeDC
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostiorID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostID] <- hostiorID
    
    
    
    ###calculate acceptance probability
    logproposalratio <- log(intervalDOWN/intervalNOW)
    #     print(paste0("logpropratio = ",logproposalratio))
    
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
    #     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  
  
  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}

.updatepathUPDOWN <- function() {
  evalq({
    ###change to proposal state
    hostioriorID <- v$nodehosts[hostiorID + 2*p$obs - 1]
    
    #for both the focal host and its infector
    v$nodetimes[hostID + 2*p$obs - 1] <- v$nodetimes[hostiorID + 2*p$obs - 1]
    v$nodetimes[hostiorID + 2*p$obs - 1] <- tinf.prop
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    v$nodeparents[nodeDP] <- 2*p$obs - 1 + hostID
    v$nodeparents[2*p$obs - 1 + hostID] <- nodeUP
    v$nodeparents[nodeUDC] <- 2*p$obs - 1 + hostiorID
    v$nodeparents[2*p$obs - 1 + hostiorID] <- nodeUC
    
    for(i in nodesinvolved) {
      if(any(.ptr(v$nodeparents, i) == 2*p$obs - 1 + hostiorID)) {
        v$nodehosts[i] <- hostiorID
      } else {
        v$nodehosts[i] <- hostID
      }
    }
    v$nodehosts[2*p$obs - 1 + hostiorID] <- hostID
    v$nodehosts[2*p$obs - 1 + hostID] <- hostioriorID
    
    
    
    ###calculate acceptance probability
    logproposalratio <- log(intervalUPDOWN/intervalNOW)
    #     print(paste0("logpropratio = ",logproposalratio))
    
  },
  .phybreakenv.prop)
  
  .propose.phybreakenv() 
  
  evalq({
    logaccprob <- logLikgen +
      logLiksam + logLikcoal - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam - .phybreakenv$logLikcoal + logproposalratio
    #     print(paste0("logaccprob = ",logaccprob))
  },.phybreakenv.prop)
  
  
  
  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    #    print("accepted")
    .accept.phybreakenv()
  }
  
  
}



}


###For updating the parameters
.update.mu <- function() {
  .prepare.phybreakenv()
  
  evalq(
    p$mu<-exp(log(p$mu) +
                rnorm(1,0,2.38*sd(log(h$si),na.rm=TRUE))),
    .phybreakenv.prop
  )


  .propose.phybreakenv()

  evalq(
    logaccprob <- logLik + logLikgen +
      logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam,
    .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    .accept.phybreakenv()
  }

  evalq(h$si <- c(h$si[-1], p$mu), .phybreakenv)
}



.update.mS <- function() {
  evalq({
    sumst <- sum(v$nodetimes[v$nodetypes == "s"] -
                        v$nodetimes[v$nodetypes == "t"])
    p$mean.sample <- p$shape.sample / rgamma(1,shape = p$shape.sample * p$obs + h$mS.sh,
                                             rate = sumst + h$mS.sh * h$mS.av/p$shape.sample)
    logLiksam <-
      .lik.sampletimes(p$shape.sample,
                       p$mean.sample,
                       v$nodetimes,
                       v$nodetypes)
    rm(sumst)
  },.phybreakenv)


}

.update.mG <- function() {
  evalq({
    sumgt <- sum(v$nodetimes[v$nodetypes == "t" & v$nodehosts != 0] -
                   v$nodetimes[v$nodetypes == "t"][v$nodehosts[v$nodetypes == "t"]])
    p$mean.gen <- p$shape.gen / rgamma(1,shape = p$shape.gen * (p$obs-1) + h$mG.sh,
                                          rate = sumgt + h$mG.sh * h$mG.av/p$shape.gen)
    logLikgen <-
      .lik.gentimes(p$obs, p$shape.gen,
                    p$mean.gen, v$nodetimes,
                    v$nodehosts, v$nodetypes
      )
    rm(sumgt)
  },
  .phybreakenv)

}

.update.wh <- function() {
  
  logLikcurrent <- 
    .lik.coaltimes(.phybreakenv$p$obs, .phybreakenv$p$wh.model, .phybreakenv$p$wh.rate0, .phybreakenv$p$wh.lambda,
                   .phybreakenv$p$wh.slope, .phybreakenv$v$nodetimes, .phybreakenv$v$nodehosts, .phybreakenv$v$nodetypes)
  if(logLikcurrent > -Inf) {
    #  print(paste0("slope = ",.phybreakenv$p$wh.slope,"; lLc = ",logLikcurrent))
    if(.phybreakenv$p$wh.model == 3) {
      proposerate0 <- exp(log(.phybreakenv$p$wh.rate0) +
                            rnorm(1,0,2.38*sd(log(.phybreakenv$h$si.wh),na.rm=TRUE)))
      proposeslope <- .phybreakenv$p$wh.slope
    } else {
      proposerate0 <- .phybreakenv$p$wh.rate0
      proposeslope <- exp(log(.phybreakenv$p$wh.slope) +
                            rnorm(1,0,2.38*sd(log(.phybreakenv$h$si.wh),na.rm=TRUE)))
    }
    logLikpropose <- 
      .lik.coaltimes(.phybreakenv$p$obs, .phybreakenv$p$wh.model, proposerate0, .phybreakenv$p$wh.lambda,
                     proposeslope, .phybreakenv$v$nodetimes, .phybreakenv$v$nodehosts, .phybreakenv$v$nodetypes)
    #  print(paste0("lLp = ",logLikpropose))
    logaccprob <- logLikpropose - logLikcurrent + 
      log(proposerate0) + log(proposeslope) -
      log(.phybreakenv$p$wh.rate0) - log(.phybreakenv$p$wh.slope) +
      dgamma(proposerate0, shape = .phybreakenv$h$wh.sh, 
             scale = .phybreakenv$h$wh.av/.phybreakenv$h$wh.sh, log = TRUE) - 
      dgamma(.phybreakenv$p$wh.rate0, shape = .phybreakenv$h$wh.sh, 
             scale = .phybreakenv$h$wh.av/.phybreakenv$h$wh.sh, log = TRUE) +
      dgamma(proposeslope, shape = .phybreakenv$h$wh.sh, 
             scale = .phybreakenv$h$wh.av/.phybreakenv$h$wh.sh, log = TRUE) - 
      dgamma(.phybreakenv$p$wh.slope, shape = .phybreakenv$h$wh.sh, 
             scale = .phybreakenv$h$wh.av/.phybreakenv$h$wh.sh, log = TRUE)
    #   print(c(.phybreakenv$p$wh.rate0,proposerate0,
    #           .phybreakenv$logLikcoal,logLikpropose,logaccprob))
    if(runif(1,0,1) < exp(logaccprob)) {
      .phybreakenv$p$wh.rate0 <- proposerate0
      .phybreakenv$p$wh.slope <- proposeslope
      .phybreakenv$logLikcoal <- logLikpropose
      #   print("accepted")
    } else {
      .phybreakenv$logLikcoal <- logLikcurrent
    }
    
    evalq(h$si.wh <- c(h$si.wh[-1], if(p$wh.model == 3) {p$wh.rate0} else {p$wh.slope}), .phybreakenv)
    
  }
  
}



