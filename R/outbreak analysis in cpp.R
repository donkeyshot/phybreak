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

#The environments are only used during MCMC-updating.

.phybreakenv <- new.env()
.phybreakenv.prop <- new.env()


.build.phybreakenv <- function(phybreak.obj) {
  .phybreakenv.prop$d <- phybreak.obj$d
  .phybreakenv.prop$v <- phybreak.obj$v
  .phybreakenv.prop$p <- phybreak.obj$p
  .phybreakenv.prop$h <- phybreak.obj$h
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
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "a")
  .phybreakenv.prop$likarray[cbind(2,
                             rep(1:length(.phybreakenv.prop$d$SNPfr),
                                 .phybreakenv.prop$p$obs),
                             rep(1:.phybreakenv.prop$p$obs,
                                 each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "c")
  .phybreakenv.prop$likarray[cbind(3,
                             rep(1:length(.phybreakenv.prop$d$SNPfr),
                                 .phybreakenv.prop$p$obs),
                             rep(1:.phybreakenv.prop$p$obs,
                                 each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "g")
  .phybreakenv.prop$likarray[cbind(4,
                             rep(1:length(.phybreakenv.prop$d$SNPfr),
                                 .phybreakenv.prop$p$obs),
                             rep(1:.phybreakenv.prop$p$obs,
                                 each = length(.phybreakenv.prop$d$SNPfr))
  )
  ] <- 1 * t(.phybreakenv.prop$d$SNP == "t")


  .likseqenv(.phybreakenv.prop,
             (.phybreakenv.prop$p$obs + 1):(2*.phybreakenv.prop$p$obs - 1),
            1:.phybreakenv.prop$p$obs)

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


}



logLik.phybreak <- function(phybreak.obj, seq = TRUE, t.gen = FALSE,
                            t.sample = FALSE) {
  if(seq) {
    with(phybreak.obj, .likseq(t(d$SNP), d$SNPfr,
         v$nodeparents, v$nodetimes, p$mu,p$obs))
  } else 0 + if(t.gen) {
    with(phybreak.obj, .lik.gentimes(p$obs, p$shape.gen, p$mean.gen,
         v$nodetimes, v$nodehosts, v$nodetypes))
  } else 0 + if(t.sample) {
    with(phybreak.obj, .lik.sampletimes(p$shape.sample, p$mean.sample,
         v$nodetimes, v$nodetypes))
  } else 0
}


###For updating the global outbreak variables
burnin.phybreak <- function(phybreak.object, nburnin) {
#   parked.samples <- phybreak.object$s
#   res <- phybreak.object
#   res$s <- c()
  .build.phybreakenv(phybreak.object)


  for(rep in 1:nburnin) {
    for(i in sample(phybreak.object$p$obs)) {
      .updatehost(i)
    }
    if(phybreak.object$h$est.mG) .update.mG()
    if(phybreak.object$h$est.mS) .update.mS()
    .update.mu()
  }

  res <- .destroy.phybreakenv(phybreak.object$s)



  return(res)
}

sample.phybreak <- function(phybreak.object, nsample, thin) {
  s.post <- list(
    nodetimes = with(phybreak.object,cbind(s$nodetimes, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodehosts = with(phybreak.object,cbind(s$nodehosts, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
    nodeparents = with(phybreak.object,cbind(s$nodeparents, matrix(NA, nrow=3*p$obs - 1, ncol=nsample))),
    mu = c(phybreak.object$s$mu, rep(NA, nsample)),
    mG = c(phybreak.object$s$mG, rep(NA, nsample)),
    mS = c(phybreak.object$s$mS, rep(NA, nsample)),
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
      .update.mu()
    }
    s.post$nodetimes[, sa] <- tail(.phybreakenv$v$nodetimes, -phybreak.object$p$obs)
    s.post$nodehosts[, sa] <- tail(.phybreakenv$v$nodehosts, -phybreak.object$p$obs)
    s.post$nodeparents[, sa] <- .phybreakenv$v$nodeparents
    s.post$mu[sa] <- .phybreakenv$p$mu
    s.post$mG[sa] <- .phybreakenv$p$mean.gen
    s.post$mS[sa] <- .phybreakenv$p$mean.sample
    s.post$logLik[sa] <- .phybreakenv$logLik + .phybreakenv$logLiksam + .phybreakenv$logLikgen
  }

  res <- .destroy.phybreakenv(s.post)


  return(res)

}


.updatehost <- function(hostID) {
  .prepare.phybreakenv()
  .phybreakenv.prop$hostID <- hostID
  evalq({
    #    tinf.prop <- v$nodetimes[hostID] - rgamma(1,shape=p$shape.sample,scale=p$mean.sample/p$shape.sample)
    tinf.prop <- v$nodetimes[hostID] - rgamma(1,shape=2,scale=p$mean.sample/2)

    ##going down the decision tree
    if(tinf.prop < suppressWarnings(min(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID]))) {
      #Q1=Y : tinf.prop before first infectee
      if(v$nodehosts[hostID + 2*p$obs - 1] == 0) {
        #Q12=YY : tinf.prop before first infectee & hostID is index
        .updatepathA()
      } else {
        #Q12=YN : tinf.prop before first infectee & hostID is not index
        if(tinf.prop < v$nodetimes[v$nodehosts == 0]) {
          #Q123=YNY : hostID is not index & tinf.prop is before infection of index
          iees <- which(v$nodehosts == (which(v$nodehosts == 0) - 2*p$obs + 1) & v$nodetypes == "t")
          if(hostID == iees[order(v$nodetimes[iees])][1] - 2*p$obs + 1) {
            #Q1234=YNYY : tinf.prop before infection of index & hostID is index's first infectee
            .updatepathB()
          } #Q1234=YNYN : tinf.prop before infection of index & hostID is not index nor index's first infectee
            ### DO NOTHING
        } else {
          #Q123=YNN : tinf.prop before first infectee and after infection of index
          .updatepathC()
        }
      }
    } else
    {
      #Q1=N : tinf.prop after first infectee
      if(v$nodehosts[hostID + 2*p$obs - 1] == 0) {
        #Q12=NY : tinf.prop after first infectee & hostID is index
        if(tinf.prop <= sort(tail(v$nodetimes,p$obs))[3]) {
          #Q123=NYY : hostID is index & tinf.prop is after first infectee but before third infection in outbreak
          ### DO NOTHING
          .updatepathD()
        } #Q123=NYN : hostID is index & tinf.prop is after third infection in outbreak
      } else {
        #Q12=NN : tinf.prop after first infectee & hostID is not index
        .updatepathE()
      }
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
                       p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times

    v$nodeparents[v$nodehosts == hostID] <-     #change the parent nodes of all nodes in hostID...
      .sampletopology(which(v$nodehosts == hostID), v$nodetimes[v$nodehosts == hostID], v$nodetypes[v$nodehosts == hostID],
                      hostID + 2*p$obs - 1, p$wh.model)
    #...to a correct topology, randomized where possible



    ###calculate acceptance probability
    logproposalratio <-
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = 2, scale = p$mean.sample/2,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = 2, scale = p$mean.sample/2, log=TRUE)

  },
  .phybreakenv.prop)






  .propose.phybreakenv()

  evalq(
  logaccprob <- logLik + logLikgen +
    logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
    .phybreakenv$logLiksam + logproposalratio,
  .phybreakenv.prop)



  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    .accept.phybreakenv()
  }

}

.updatepathB <- function() {

  evalq({
    ###change to proposal state
    index.current.ID <- which(v$nodehosts == 0) - 2*p$obs + 1


    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop

    v$nodehosts[hostID + 2*p$obs - 1] <- 0
    v$nodeparents[hostID + 2*p$obs - 1] <- 0
    v$nodehosts[index.current.ID + 2*p$obs - 1] <- hostID
    v$nodehosts[v$nodehosts == index.current.ID & v$nodetypes == "c"][1] <- hostID

    for(ID in c(hostID,index.current.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible

    }


    ###calculate acceptance probability
    logproposalratio <-
      dgamma(.phybreakenv$v$nodetimes[hostID] -
               .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
             shape = 2, scale = p$mean.sample/2,log=TRUE) -
      dgamma(v$nodetimes[hostID] -
               v$nodetimes[hostID + 2*p$obs - 1],
             shape = 2, scale = p$mean.sample/2, log=TRUE)

  },.phybreakenv.prop)




  .propose.phybreakenv()

  evalq(
    logaccprob <- logLik + logLikgen +
      logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam + logproposalratio,
    .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
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
                         p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times

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
           shape = 2, scale = p$mean.sample/2,log=TRUE) -
    dgamma(v$nodetimes[hostID] -
             v$nodetimes[hostID + 2*p$obs - 1],
           shape = 2, scale = p$mean.sample/2, log=TRUE)

  },.phybreakenv.prop)




  .propose.phybreakenv()

  evalq(
  logaccprob <- logLik + logLikgen +
    logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
    .phybreakenv$logLiksam + logproposalratio,
  .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    .accept.phybreakenv()
  }

}

.updatepathD <- function() {


  evalq({
    ###change to proposal state

    iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
    infectee.current.ID <-
      iees.nodeIDs[v$nodetimes[iees.nodeIDs] ==
                     min(v$nodetimes[iees.nodeIDs])] - 2*p$obs + 1
    v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop

    v$nodehosts[infectee.current.ID + 2*p$obs - 1] <- 0
    v$nodeparents[infectee.current.ID + 2*p$obs - 1] <- 0
    v$nodehosts[hostID + 2*p$obs - 1] <- infectee.current.ID
    v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infectee.current.ID

    for(ID in c(hostID,infectee.current.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
        v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
                         p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible

    }

  ###calculate acceptance probability
  logproposalratio <-
    dgamma(.phybreakenv$v$nodetimes[hostID] -
             .phybreakenv$v$nodetimes[hostID + 2*p$obs - 1],
           shape = 2, scale = p$mean.sample/2,log=TRUE) -
    dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
           shape = 2, scale = p$mean.sample/2, log=TRUE)

  },.phybreakenv.prop)





  .propose.phybreakenv()

  evalq(
  logaccprob <- logLik + logLikgen +
    logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
    .phybreakenv$logLiksam + logproposalratio,
  .phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
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
                         p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times

      v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
        .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
                        ID + 2*p$obs - 1, p$wh.model)
      #...to a correct topology, randomized where possible

    }

  ###calculate acceptance probability

  logproposalratio <-
    pgamma(v$nodetimes[infectee.current.ID] - v$nodetimes[infectee.current.ID + 2*p$obs - 1],
           shape = 2, scale = p$mean.sample/2, log.p = TRUE) -
    pgamma(v$nodetimes[hostID] - v$nodetimes[infector.current.ID + 2*p$obs - 1],
           shape = 2, scale = p$mean.sample/2, log.p = TRUE)

  },.phybreakenv.prop)




  .propose.phybreakenv()

  evalq({
    logaccprob <- logLik + logLikgen +
      logLiksam - .phybreakenv$logLik - .phybreakenv$logLikgen -
      .phybreakenv$logLiksam + logproposalratio
  },.phybreakenv.prop)


  if(runif(1,0,1) < exp(.phybreakenv.prop$logaccprob)) {
    .accept.phybreakenv()
  }


}

}

###For updating the parameters
.update.mu <- function() {

  evalq(
    p$mu<-exp(log(p$mu) +
                rnorm(1,0,sd(log(h$si),na.rm=TRUE))),
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
                                             rate = sumst + h$mS.sc/p$shape.sample)
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
                                          rate = sumgt + h$mG.sc/p$shape.gen)
    logLikgen <-
      .lik.gentimes(p$obs, p$shape.gen,
                    p$mean.gen, v$nodetimes,
                    v$nodehosts, v$nodetypes
      )
    rm(sumgt)
  },
  .phybreakenv)

}




