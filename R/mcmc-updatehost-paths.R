### all helper functions exclusively involved in updating the tree, ### i.e. proposing a new tree and accepting or rejecting
### the proposal ### fixed parameters
tinf.prop.shape.mult <- 2/3  #shape for proposing infection time is shape.sample * tinf.prop.shape.mult

### updating both transmission and phylotree: proposing an infection time and following the decision tree called from:
### burnin.phybreak sample.phybreak calling: .prepare.pbe .updatepathA, .updatepathB, .updatepathC1, .updatepathC2,
### .updatepathD, .updatepathE, .updatepathF1, .updatepathF2
.updatehost <- function(hostID) {
    ### create an up-to-date proposal-environment with hostID as focal host
    .prepare.pbe()
    .copy2pbe1("hostID", environment())
  
    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    p <- .pbe1$p
    v <- .pbe1$v
    
    ### propose the new infection time
    tinf.prop <- v$nodetimes[hostID] - 
      rgamma(1, shape = tinf.prop.shape.mult * p$shape.sample, scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample))
    .copy2pbe1("tinf.prop", le)

    ### going down the decision tree
    if (v$nodehosts[v$nodetypes == "t"][hostID] == 0) {
      # Y (hostID is index case)
      if (tinf.prop < min(c(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID], Inf))) {
        # YY (... & tinf.prop before hostID's first transmission node)
        .updatepathA()
      } else {
        # YN (... & tinf.prop after hostID's first transmission node)
        if (tinf.prop < sort(c(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID], Inf))[2]) {
          # YNY (... & tinf.prop before hostID's second transmission node)
          .updatepathB()
        } else {
          # YNN (... & tinf.prop after hostID's second transmission node)
          if (runif(1) < 0.5) 
            .updatepathC(TRUE) else .updatepathC(FALSE)
        }
      }
    } else {
      # N (hostID is not index case)
      if (tinf.prop < v$nodetimes[v$nodehosts == 0]) {
        # NY (... & tinf.prop before infection of index case)
        .updatepathD()
      } else {
        # NN (... & tinf.prop after infection of index case)
        if (tinf.prop < min(c(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID], Inf))) {
          # NNY (... & tinf.prop before hostID's first transmission node)
          .updatepathE()
        } else {
          # NNN (... & tinf.prop after hostID's first transmission node)
          if (runif(1) < 0.5) 
            .updatepathF(TRUE) else .updatepathF(FALSE)
        }
      }
    }
}


### updating transmission tree but keeping phylotree: proposing an infection time and following the decision tree called from:
### burnin.phybreak sample.phybreak calling: .prepare.pbe .ptr ##C++ .updatepathG, .updatepathH, .updatepathI,
### .updatepathJ
.updatehost.keepphylo <- function(hostID) {
    ### create an up-to-date proposal-environment with hostID as focal host
    .prepare.pbe()
    .copy2pbe1("hostID", environment())
  
    ### making variables and parameters available within the function
    le <- environment()
    p <- .pbe1$p
    v <- .pbe1$v
  
    ### propose the new infection time
    tinf.prop <- v$nodetimes[hostID] - 
      rgamma(1, shape = tinf.prop.shape.mult * p$shape.sample, scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample))
    .copy2pbe1("tinf.prop", le)

    ### identify the focal host's infector
    hostiorID <- v$nodehosts[2 * p$obs - 1 + hostID]
    .copy2pbe1("hostiorID", le)
    
    ### going down the decision tree
    if (hostiorID == 0) {
      # Y (hostID is index)
      if (tinf.prop < v$nodetimes[which(v$nodeparents == 2 * p$obs - 1 + hostID)]) {
        # YY (... & tinf.prop before root of phylotree)
        .updatepathG()
      } else {
        # YN (... & tinf.prop after root of phylotree): reject
      }
    } else {
      # N (hostID is not index)
      nodemrca <- .ptr(v$nodeparents, hostID)[.ptr(v$nodeparents, hostID) %in% .ptr(v$nodeparents, hostiorID)][1]
      timemrca <- v$nodetimes[nodemrca]
      .copy2pbe1("nodemrca", le)
      .copy2pbe1("timemrca", le)
      if (tinf.prop > timemrca) {
        # NY (... & tinf.prop after MRCA of hostID and hostiorID)
        .updatepathH()
      } else {
        # NN (... & tinf.prop before MRCA of hostID and hostiorID)
        tinf2.prop <- v$nodetimes[hostiorID] - 
          rgamma(1, shape = tinf.prop.shape.mult * p$shape.sample, scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample))
        .copy2pbe1("tinf2.prop", le)
        if (tinf2.prop > timemrca) {
          # NNY (... & tinf2.prop after MRCA of hostID and hostiorID)
          hostioriorID <- v$nodehosts[2 * p$obs - 1 + hostiorID]
          .copy2pbe1("hostioriorID", le)
          if (hostioriorID == 0) {
            # NNYY (... & hostiorID is index)
            tinf.prop <- v$nodetimes[hostiorID + 2 * p$obs - 1]
            .copy2pbe1("tinf.prop", le)
            .updatepathI()
          } else {
            # NNYN (... & hostiorID is not index)
            nodemrcaior <- .ptr(v$nodeparents, hostID)[.ptr(v$nodeparents, hostID) %in% .ptr(v$nodeparents, hostioriorID)][1]
            timemrcaior <- v$nodetimes[nodemrcaior]
            .copy2pbe1("nodemrcaior", le)
            .copy2pbe1("timemrcaior", le)
            if (tinf.prop > timemrcaior) {
              # NNYNY (... & tinf.prop after MRCA of hostID and hostioriorID)
              .updatepathJ()
            } else {
              # NNYNN (... & tinf.prop before MRCA of hostID and hostioriorID): reject
            }
          }
        } else {
          # NNN (... & tinf2.prop before MRCA of hostID and hostiorID): reject
        }
      }
    }
}


{
  ### update if hostID is index and tinf.prop is before the first secondary case called from: .updatehost calling:
  ### .samplecoaltimes .sampletopology .propose.pbe .accept.pbe
  .updatepathA <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
    hostID <- .pbe1$hostID
    tinf.prop <- .pbe1$tinf.prop
    
    ### change to proposal state
    
    # bookkeeping: change infection time
    v$nodetimes[v$nodetypes == "t"][hostID] <- tinf.prop
    
    # phylotree in hostID
    v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-
      v$nodetimes[v$nodetypes == "t"][hostID] +
      .samplecoaltimes(v$nodetimes[v$nodehosts == hostID & v$nodetypes != "c"] - v$nodetimes[v$nodetypes == "t"][hostID],
                       p$wh.model, p$wh.slope)
    v$nodeparents[v$nodehosts == hostID] <-
      .sampletopology(which(v$nodehosts == hostID),
                      v$nodetimes[v$nodehosts == hostID],
                      v$nodetypes[v$nodehosts == hostID],
                      hostID + 2 * d$nsamples - 1, p$wh.model)
    
    ### update proposal environment
    .copy2pbe1("v", le)
    
    ### calculate proposal ratio
    logproposalratio <- dgamma(.pbe0$v$nodetimes[.pbe1$hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"][.pbe1$hostID],
                               shape = tinf.prop.shape.mult * p$shape.sample,
                               scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) -
      dgamma(v$nodetimes[hostID] - v$nodetimes[v$nodetypes == "t"][hostID],
             shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)
    
    
    ### calculate likelihood
    .propose.pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikseq + .pbe1$logLikgen + .pbe1$logLiksam - .pbe0$logLikseq - .pbe0$logLikgen - .pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      .accept.pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is index and tinf.prop is after the first secondary case, but before the second secondary case called
  ### from: .updatehost calling: .samplecoaltimes .sampletopology .propose.pbe .accept.pbe
  .updatepathB <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
    hostID <- .pbe1$hostID
    tinf.prop <- .pbe1$tinf.prop
    
    ### change to proposal state
    
    # identify the current first infectee of hostID
    iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
    infectee.current.ID <- iees.nodeIDs[v$nodetimes[iees.nodeIDs] == min(v$nodetimes[iees.nodeIDs])] - 2 * d$nsamples + 1
    
    # bookkeeping: change infection time
    v$nodetimes[v$nodetypes == "t"][hostID] <- tinf.prop
    
    # propose infector for hostID
    dens.infectorproposal <- dgamma(tinf.prop - v$nodetimes[v$nodetypes == "t"],
                                    shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) +
      (tinf.prop - v$nodetimes[v$nodetypes == "t"] > 0)/.pbe0$h$dist[hostID, ]
    dens.infectorproposal[hostID] <- 0
    infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
    
    # prepare for phylotree proposals by moving the coalescent node from hostID to the new infector
    ## step 1: remove moving coalescent node from hostID, and place it before hostID
    movingcoalnode <- v$nodeparents[v$nodetypes == "t"][infectee.current.ID] #identify moving coalescent node
    otherupnode <- setdiff(which(v$nodeparents == movingcoalnode), 2 * d$nsamples + infectee.current.ID - 1) #identify second childnode of moving node
    v$nodeparents[v$nodetypes == "t"][infectee.current.ID] <- 0 #place rootnode before new index case
    v$nodeparents[v$nodetypes == "t"][hostID] <- movingcoalnode #place moving node before hostID
    v$nodeparents[otherupnode] <- v$nodeparents[movingcoalnode] #remove moving node from hostID
    ## step 2: sample new coalescent time and edge within which moving node should be placed
    v$nodetimes[movingcoalnode] <- v$nodetimes[v$nodetypes == "t"][infector.proposed.ID] +
      .sampleextracoaltime(v$nodetimes[v$nodehosts == infector.proposed.ID] - v$nodetimes[v$nodetypes == "t"][infector.proposed.ID],
                           v$nodetypes[v$nodehosts == infector.proposed.ID],
                           v$nodetimes[v$nodetypes == "t"][hostID] - v$nodetimes[v$nodetypes == "t"][infector.proposed.ID],
                           p$wh.model, p$wh.slope) #sample new coalescent time
    otherupnode <- .sampleextraupnode(which(v$nodehosts == infector.proposed.ID),
                                      v$nodeparents[v$nodehosts == infector.proposed.ID],
                                      v$nodetimes[v$nodehosts == infector.proposed.ID],
                                      v$nodetimes[movingcoalnode])  #sample edge at which moving node should be attached
    ## step 3: place moving coalescent node inside proposed infector
    v$nodeparents[movingcoalnode] <- v$nodeparents[otherupnode] #attach moving node to minitree in proposed infector
    v$nodeparents[otherupnode] <- movingcoalnode #attacb minitree to moving node
    v$nodehosts[movingcoalnode] <- infector.proposed.ID #place moving node inside proposed infector
    
    # bookkeeping: change the transmission tree topology
    v$nodehosts[v$nodetypes == "t"][infectee.current.ID] <- 0
    v$nodehosts[v$nodetypes == "t"][hostID] <- infector.proposed.ID
    
    # propose phylotrees for hostID and the new infector
    for (ID in c(hostID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[v$nodetypes == "t"][ID] +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[v$nodetypes == "t"][ID],
                         p$wh.model, p$wh.slope)
      v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID],
                                                          v$nodetypes[v$nodehosts == ID], ID + 2 * d$nsamples - 1, p$wh.model)
    }
    
    
    # add coalescent node to the new infector
    ## in: cbind(nodeparents[nodehosts == infector.proposed.ID], which(nodehosts == infector.proposed.ID)) #edges
    ## in: nodetimes[nodehosts == infector.proposed.ID] - nodetimes[nodetypes=="t][infector.proposed.ID] #times of edgeends
    
    ### update proposal environment
    .copy2pbe1("v", le)
    
    ### calculate proposal ratio
    logproposalratio <- log(sum(dens.infectorproposal)/(dens.infectorproposal[infector.proposed.ID])) +
      dgamma(.pbe0$v$nodetimes[.pbe1$hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"][.pbe1$hostID],
             shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) -
      dgamma(v$nodetimes[hostID] - v$nodetimes[v$nodetypes == "t"][hostID],
             shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)
    
    ### calculate likelihood
    .propose.pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikseq + .pbe1$logLikgen + .pbe1$logLiksam - .pbe0$logLikseq - .pbe0$logLikgen - .pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      .accept.pbe("phylotrans")
    }
    
  }
  
  
  ### update if hostID is index and tinf.prop is after the second secondary case called from: .updatehost calling:
  ### .samplecoaltimes .sampletopology .propose.pbe .accept.pbe
  .updatepathC <- function(exchange) {
    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
    hostID <- .pbe1$hostID
    tinf.prop <- .pbe1$tinf.prop
    
    ### change to proposal state
    
    # identify the current first infectee of hostID
    iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
    infectee.current.ID <- iees.nodeIDs[v$nodetimes[iees.nodeIDs] == min(v$nodetimes[iees.nodeIDs])] - 2 * d$nsamples + 1

    # identify the other infectees of hostID and those of infectee.current.ID
    infecteenodes.primary <- iees.nodeIDs[iees.nodeIDs != (infectee.current.ID + 2 * d$nsamples - 1)]
    infecteenodes.secondary <- which(v$nodehosts == infectee.current.ID & v$nodetypes == "t")

    # remember current infector and infection time of hostID
    infector.current.ID <- v$nodehosts[v$nodetypes == "t"][hostID]
    inftime.current.hostID <- v$nodetimes[v$nodetypes == "t"][hostID]

    # bookkeeping: change infection times
    v$nodetimes[v$nodetypes == "t"][hostID] <- v$nodetimes[v$nodetypes == "t"][infectee.current.ID]
    v$nodetimes[v$nodetypes == "t"][infectee.current.ID] <- inftime.current.hostID
    
    # bookkeeping: change the transmission tree topology
    v$nodehosts[v$nodetypes == "t"][hostID] <- infectee.current.ID
    v$nodehosts[v$nodetypes == "t"][infectee.current.ID] <- infector.current.ID
    v$nodeparents[v$nodetypes == "t"][infectee.current.ID] <- v$nodeparents[v$nodetypes == "t"][hostID]
    if (exchange) {
      v$nodehosts[infecteenodes.primary] <- infectee.current.ID
      v$nodehosts[infecteenodes.secondary] <- hostID
    }

    # prepare for phylotree proposals by moving coalescent nodes between hostID and infectee.current.ID
    if (exchange) {
      coalnodes.bothhosts <- which((v$nodehosts %in% c(hostID, infectee.current.ID)) & v$nodetypes == "c")
      nrcoalnodes.infecteeID <- sum(v$nodehosts == infectee.current.ID & v$nodetypes != "c") - 1
      v$nodehosts[head(coalnodes.bothhosts, nrcoalnodes.infecteeID)] <- infectee.current.ID
      v$nodehosts[tail(coalnodes.bothhosts, -nrcoalnodes.infecteeID)] <- hostID
    } else {
      v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infectee.current.ID
    }

    # propose phylotrees for hostID and the new infector
    for (ID in c(hostID, infectee.current.ID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[v$nodetypes == "t"][ID] +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[v$nodetypes == "t"][ID],
                         p$wh.model, p$wh.slope)
      v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID],
                                                          v$nodetypes[v$nodehosts == ID], ID + 2 * d$nsamples - 1, p$wh.model)
    }
    
    ### update proposal environment
    .copy2pbe1("v", le)
    
    ### calculate proposal ratio
    if (infector.current.ID != 0) {
      ## .updatepathF: probability to sample after FIRST infectee
      logproposalratio <- pgamma(v$nodetimes[infectee.current.ID] - v$nodetimes[v$nodetypes == "t"][hostID],
                                 shape = tinf.prop.shape.mult * p$shape.sample,
                                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log.p = TRUE) -
        pgamma(v$nodetimes[hostID] - v$nodetimes[v$nodetypes == "t"][hostID],
               shape = tinf.prop.shape.mult * p$shape.sample,
               scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log.p = TRUE)
    } else {
      ## .updatepathC: probability to sample after SECOND infectee
      
      # second infectee not necessarily the same for reversal proposal, so identify these intervals for hostID and new index
      sampleinterval.current.ID <- v$nodetimes[hostID] - sort(.pbe0$v$nodetimes[iees.nodeIDs])[2]
      iees.newindex <- which(v$nodehosts == infectee.current.ID & v$nodetypes == "t")
      sampleinterval.newindex <- v$nodetimes[infectee.current.ID] - sort(c(v$nodetimes[iees.newindex], Inf))[2]
      
      logproposalratio <- pgamma(sampleinterval.newindex, shape = tinf.prop.shape.mult * p$shape.sample,
                                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log.p = TRUE) -
        pgamma(sampleinterval.current.ID, shape = tinf.prop.shape.mult * p$shape.sample,
               scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log.p = TRUE)
    }
    
    
    ### calculate likelihood
    .propose.pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikseq + .pbe1$logLikgen + .pbe1$logLiksam - .pbe0$logLikseq - .pbe0$logLikgen - .pbe0$logLiksam +
      logproposalratio
    
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      .accept.pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is before infection of the index case called from: .updatehost calling:
  ### .samplecoaltimes .sampletopology .propose.pbe .accept.pbe
  .updatepathD <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
    hostID <- .pbe1$hostID
    tinf.prop <- .pbe1$tinf.prop
    
    
    ### change to proposal state
    
    # identify the current index and current infector
    index.current.ID <- which(v$nodehosts == 0) - 2 * d$nsamples + 1
    infector.current.ID <- v$nodehosts[v$nodetypes == "t"][hostID]
    
    # bookkeeping: change infection time
    v$nodetimes[v$nodetypes == "t"][hostID] <- tinf.prop
    
    # prepare for phylotree proposals by moving coalescent nodes from current infector to hostID
    ## step 1: remove moving coalescent node from current infector, and place it before current index
    movingcoalnode <- v$nodeparents[v$nodetypes == "t"][hostID] #identify moving coalescent node
    otherupnode <- setdiff(which(v$nodeparents == movingcoalnode), 2 * d$nsamples + hostID - 1) #identify second childnode of moving node
    v$nodeparents[otherupnode] <- v$nodeparents[movingcoalnode] #remove coalescent node from current infector's minitree
    v$nodeparents[v$nodetypes == "t"][hostID] <- 0 #place rootnode before new index case
    v$nodeparents[v$nodetypes == "t"][index.current.ID] <- movingcoalnode #place moving node before current index
    ## step 2: place moving coalescent node inside proposed index
    v$nodehosts[movingcoalnode] <- hostID #place moving node inside new index
    
    # bookkeeping: change the transmission tree topology
    v$nodehosts[v$nodetypes == "t"][hostID] <- 0
    v$nodehosts[v$nodetypes == "t"][index.current.ID] <- hostID
    
    # propose phylotrees for hostID and the new infector
    for (ID in c(hostID)) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[v$nodetypes == "t"][ID] +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[v$nodetypes == "t"][ID],
                         p$wh.model, p$wh.slope)
      v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID],
                                                          v$nodetypes[v$nodehosts == ID], ID + 2 * d$nsamples - 1, p$wh.model)
    }
    
    ### update proposal environment
    .copy2pbe1("v", le)
    
    ### calculate proposal ratio the reverse proposal includes proposing an infector
    dens.infectorcurrent <- dgamma(.pbe0$v$nodetimes[v$nodetypes == "t"][hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"],
                                   shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) +
      (.pbe0$v$nodetimes[v$nodetypes == "t"][hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"] > 0)/h$dist[hostID, ]
    dens.infectorcurrent[hostID] <- 0
    
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID]/(sum(dens.infectorcurrent))) +
      dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"][hostID],
             shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) -
      dgamma(v$nodetimes[hostID] - v$nodetimes[v$nodetypes == "t"][hostID],
             shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)
    
    
    ### calculate likelihood
    .propose.pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikseq + .pbe1$logLikgen + .pbe1$logLiksam - .pbe0$logLikseq - .pbe0$logLikgen - .pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      .accept.pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is after infection of the index case, but before the first secondary case
  ### called from: .updatehost calling: .samplecoaltimes .sampletopology .propose.pbe .accept.pbe
  .updatepathE <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- .pbe0$d
    h <- .pbe0$h
    p <- .pbe1$p
    v <- .pbe1$v
    hostID <- .pbe1$hostID
    tinf.prop <- .pbe1$tinf.prop
    
    
    ### change to proposal state
    
    # identify the current infector
    infector.current.ID <- v$nodehosts[v$nodetypes == "t"][hostID]
    
    # bookkeeping: change infection time
    v$nodetimes[v$nodetypes == "t"][hostID] <- tinf.prop
    
    # propose infector for hostID
    dens.infectorproposal <- dgamma(tinf.prop - v$nodetimes[v$nodetypes == "t"],
                                    shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) +
      (tinf.prop - v$nodetimes[v$nodetypes == "t"] > 0)/h$dist[hostID, ]
    dens.infectorproposal[hostID] <- 0
    infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
    
    # prepare for phylotree proposals by moving the coalescent node from the current infector to the proposed infector
    ## step 1: remove transmission node and moving coalescent node from current infector and from tree
    movingcoalnode <- v$nodeparents[v$nodetypes == "t"][hostID] #identify moving coalescent node
    otherupnode <- setdiff(which(v$nodeparents == movingcoalnode), 2 * d$nsamples + hostID - 1) #identify second childnode of moving node
    v$nodeparents[otherupnode] <- v$nodeparents[movingcoalnode] #remove coalescent node from current infector's minitree
    v$nodehosts[v$nodetypes == "t"][hostID] <- -1L
    v$nodehosts[movingcoalnode] <- -1L
    v$nodeparents[movingcoalnode] <- -1L
    ## step 2: sample new coalescent time and edge on which moving node should be placed
    v$nodetimes[movingcoalnode] <- v$nodetimes[v$nodetypes == "t"][infector.proposed.ID] +
      .sampleextracoaltime(v$nodetimes[v$nodehosts == infector.proposed.ID] - v$nodetimes[v$nodetypes == "t"][infector.proposed.ID],
                           v$nodetypes[v$nodehosts == infector.proposed.ID],
                           v$nodetimes[v$nodetypes == "t"][hostID] - v$nodetimes[v$nodetypes == "t"][infector.proposed.ID],
                           p$wh.model, p$wh.slope) #sample new coalescent time
    otherupnode <- .sampleextraupnode(which(v$nodehosts == infector.proposed.ID),
                                      v$nodeparents[v$nodehosts == infector.proposed.ID],
                                      v$nodetimes[v$nodehosts == infector.proposed.ID],
                                      v$nodetimes[movingcoalnode])  #sample edge at which moving node should be attached
    ## step 3: place transmission node and moving coalescent node inside proposed infector
    v$nodeparents[movingcoalnode] <- v$nodeparents[otherupnode] #attach moving node to minitree in proposed infector
    v$nodeparents[otherupnode] <- movingcoalnode #attacb minitree to moving node
    v$nodehosts[movingcoalnode] <- infector.proposed.ID #place moving node inside proposed infector
    v$nodehosts[v$nodetypes == "t"][hostID] <- infector.proposed.ID
    
    # propose phylotrees for hostID
    for (ID in unique(c(hostID))) {
      v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[v$nodetypes == "t"][ID] +
        .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[v$nodetypes == "t"][ID],
                         p$wh.model, p$wh.slope)
      v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID],
                                                          v$nodetypes[v$nodehosts == ID], ID + 2 * d$nsamples - 1, p$wh.model)
    }
    
    ### update proposal environment
    .copy2pbe1("v", le)
    
    ### calculate proposal ratio the reverse proposal includes proposing an infector
    dens.infectorcurrent <- dgamma(.pbe0$v$nodetimes[v$nodetypes == "t"][hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"],
                                   shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) +
      (.pbe0$v$nodetimes[v$nodetypes == "t"][hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"] > 0)/h$dist[hostID, ]
    dens.infectorcurrent[hostID] <- 0
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal)/
                              (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) +
      dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[v$nodetypes == "t"][hostID],
             shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) -
      dgamma(v$nodetimes[hostID] - v$nodetimes[v$nodetypes == "t"][hostID], shape = tinf.prop.shape.mult * p$shape.sample,
             scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)
    
    
    ### calculate likelihood
    .propose.pbe("phylotrans")
    
    
    ### calculate acceptance probability
    logaccprob <- .pbe1$logLikseq + .pbe1$logLikgen + .pbe1$logLiksam - .pbe0$logLikseq - .pbe0$logLikgen - .pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      .accept.pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is after the first secondary case called from: .updatehost calling:
  ### .updatepathC
  .updatepathF <- function(exchange) {
    .updatepathC(exchange)
  }
  
  
  
}



{
    
    ### update if hostID is index and tinf.prop is before the first coalescence node called from: .updatehost.keepphylo calling:
    ### .propose.pbe .accept.pbe
    .updatepathG <- function() {
        ### making variables and parameters available within the function
        le <- environment()
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        tinf.prop <- .pbe1$tinf.prop
      

        ### change to proposal state
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
        
        
        ### calculate proposal ratio
        logproposalratio <- dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[hostID + 2 * p$obs - 1], 
                                   shape = tinf.prop.shape.mult * p$shape.sample, 
                                   scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)

        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate likelihood
        .propose.pbe("trans")
        
        ### calculate acceptance probability
        logaccprob <- .pbe1$logLikgen + .pbe1$logLiksam + .pbe1$logLikcoal - .pbe0$logLikgen - .pbe0$logLiksam - .pbe0$logLikcoal + 
            logproposalratio
        
        ### accept or reject
        if (runif(1) < exp(logaccprob)) {
            .accept.pbe("trans")
        }
    }
    
    
    ### update if hostID is not index and tinf.prop is after the MRCA of hostID and infector called from: .updatehost.keepphylo
    ### calling: .ptr ##C++ .propose.pbe .accept.pbe
    .updatepathH <- function() {
        ### making variables and parameters available within the function
        le <- environment()
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        hostiorID <- .pbe1$hostiorID
        tinf.prop <- .pbe1$tinf.prop
      
        ### change to proposal state
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
        
        # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
        nodeUC <- v$nodeparents[2 * p$obs - 1 + hostID]
        nodeDC <- which(v$nodeparents == 2 * p$obs - 1 + hostID)
        v$nodeparents[nodeDC] <- nodeUC
        # second, place it between the proposed upstream and downstream nodes
        nodeUP <- .ptr(v$nodeparents, hostID)[v$nodetimes[.ptr(v$nodeparents, hostID)] < tinf.prop][1]
        nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID))
        v$nodeparents[nodeDP] <- 2 * p$obs - 1 + hostID
        v$nodeparents[2 * p$obs - 1 + hostID] <- nodeUP
        
        # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
        nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostiorID)
        for (i in nodesinvolved) {
          if (any(.ptr(v$nodeparents, i) == 2 * p$obs - 1 + hostID)) {
            v$nodehosts[i] <- hostID
          } else {
            v$nodehosts[i] <- hostiorID
          }
        }
        v$nodehosts[2 * p$obs - 1 + hostID] <- hostiorID
        
        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate proposal ratio
        logproposalratio <- dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[hostID + 2 * p$obs - 1],
                                   shape = tinf.prop.shape.mult * p$shape.sample, 
                                   scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)

        
        ### calculate likelihood
        .propose.pbe("trans")
        
        ### calculate acceptance probability
        logaccprob <- .pbe1$logLikgen + .pbe1$logLiksam + .pbe1$logLikcoal - .pbe0$logLikgen - .pbe0$logLiksam - .pbe0$logLikcoal + 
            logproposalratio
        
        ### accept or reject
        if (runif(1) < exp(logaccprob)) {
            .accept.pbe("trans")
        }
    }
    
    
    ### update if hostID is not index, tinf.prop is before the MRCA of hostID and infector, and infector is index called from:
    ### .updatehost.keepphylo calling: .ptr ##C++ .propose.pbe .accept.pbe
    .updatepathI <- function() {
        ### making variables and parameters available within the function
        le <- environment()
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        hostiorID <- .pbe1$hostiorID
        hostioriorID <- .pbe1$hostioriorID
        tinf.prop <- .pbe1$tinf.prop
        tinf2.prop <- .pbe1$tinf2.prop
        timemrca <- .pbe1$timemrca
        
        ### change to proposal state
        
        ## first the infector: move transmission node downstream
        
        # bookkeeping: change infection time
        v$nodetimes[hostiorID + 2 * p$obs - 1] <- tinf2.prop
        
        # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
        nodeUC <- v$nodeparents[2 * p$obs - 1 + hostiorID]
        nodeDC <- which(v$nodeparents == 2 * p$obs - 1 + hostiorID)
        v$nodeparents[nodeDC] <- nodeUC
        # second, place it between the proposed upstream and downstream nodes
        nodeUP <- .ptr(v$nodeparents, hostiorID)[v$nodetimes[.ptr(v$nodeparents, hostiorID)] < tinf2.prop][1]
        nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostiorID))
        v$nodeparents[nodeDP] <- 2 * p$obs - 1 + hostiorID
        v$nodeparents[2 * p$obs - 1 + hostiorID] <- nodeUP
        
        # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
        nodesinvolved <- which(v$nodehosts == hostiorID | v$nodehosts == hostioriorID)
        for (i in nodesinvolved) {
          if (any(.ptr(v$nodeparents, i) == 2 * p$obs - 1 + hostiorID)) {
            v$nodehosts[i] <- hostiorID
          } else {
            v$nodehosts[i] <- hostioriorID
          }
        }
        v$nodehosts[2 * p$obs - 1 + hostiorID] <- hostioriorID
        
        ## then hostID itself: move transmission node upstream
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
        
        # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
        nodeUC <- v$nodeparents[2 * p$obs - 1 + hostID]
        nodeDC <- which(v$nodeparents == 2 * p$obs - 1 + hostID)
        v$nodeparents[nodeDC] <- nodeUC
        # second, place it between the proposed upstream and downstream nodes
        nodeUP <- c(.ptr(v$nodeparents, hostID)[v$nodetimes[.ptr(v$nodeparents, hostID)] < tinf.prop], 0)[1]
        nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID))
        v$nodeparents[nodeDP] <- 2 * p$obs - 1 + hostID
        v$nodeparents[2 * p$obs - 1 + hostID] <- nodeUP
        
        # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
        nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostioriorID)
        for (i in nodesinvolved) {
          if (any(.ptr(v$nodeparents, i) == 2 * p$obs - 1 + hostID)) {
            v$nodehosts[i] <- hostID
          } else {
            v$nodehosts[i] <- hostioriorID
          }
        }
        v$nodehosts[2 * p$obs - 1 + hostID] <- hostioriorID
        
        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate proposal ratio
        logproposalratio <- dgamma(.pbe0$v$nodetimes[.pbe1$hostID] - .pbe0$v$nodetimes[.pbe1$hostID + 2 * p$obs - 1],
                                   shape = tinf.prop.shape.mult * p$shape.sample,
                                   scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          log(1 - pgamma(v$nodetimes[hostID] - timemrca,
                         shape = tinf.prop.shape.mult * p$shape.sample, 
                         scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample))) + 
          log(1 - pgamma(.pbe0$v$nodetimes[.pbe1$hostiorID] - timemrca, 
                         shape = tinf.prop.shape.mult * p$shape.sample, 
                         scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample))) - 
          dgamma(v$nodetimes[hostiorID] - v$nodetimes[hostiorID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)

      
        
        ### calculate likelihood
        .propose.pbe("trans")
        
        ### calculate acceptance probability
        logaccprob <- .pbe1$logLikgen + .pbe1$logLiksam + .pbe1$logLikcoal - .pbe0$logLikgen - .pbe0$logLiksam - .pbe0$logLikcoal + 
            logproposalratio
        
        ### accept or reject
        if (runif(1) < exp(logaccprob)) {
            .accept.pbe("trans")
        }
    }
    
    
    ### update if hostID is not index, tinf.prop is before the MRCA of hostID and infector, infector is not index, and tinf.prop is
    ### after the MRCA of hostID and infector's infector called from: .updatehost.keepphylo calling: .ptr ##C++
    ### .propose.pbe .accept.pbe
    .updatepathJ <- function() {
        ### making variables and parameters available within the function
        le <- environment()
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        hostiorID <- .pbe1$hostiorID
        hostioriorID <- .pbe1$hostioriorID
        tinf.prop <- .pbe1$tinf.prop
        tinf2.prop <- .pbe1$tinf2.prop

        ### change to proposal state
        
        ## first the infector: move transmission node downstream
        
        # bookkeeping: change infection time
        v$nodetimes[hostiorID + 2 * p$obs - 1] <- tinf2.prop
        
        # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
        nodeUC <- v$nodeparents[2 * p$obs - 1 + hostiorID]
        nodeDC <- which(v$nodeparents == 2 * p$obs - 1 + hostiorID)
        v$nodeparents[nodeDC] <- nodeUC
        # second, place it between the proposed upstream and downstream nodes
        nodeUP <- .ptr(v$nodeparents, hostiorID)[v$nodetimes[.ptr(v$nodeparents, hostiorID)] < tinf2.prop][1]
        nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostiorID))
        v$nodeparents[nodeDP] <- 2 * p$obs - 1 + hostiorID
        v$nodeparents[2 * p$obs - 1 + hostiorID] <- nodeUP
        
        # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
        nodesinvolved <- which(v$nodehosts == hostiorID | v$nodehosts == hostioriorID)
        for (i in nodesinvolved) {
          if (any(.ptr(v$nodeparents, i) == 2 * p$obs - 1 + hostiorID)) {
            v$nodehosts[i] <- hostiorID
          } else {
            v$nodehosts[i] <- hostioriorID
          }
        }
        v$nodehosts[2 * p$obs - 1 + hostiorID] <- hostioriorID
        
        ## then hostID itself: move transmission node upstream
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
        
        # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
        nodeUC <- v$nodeparents[2 * p$obs - 1 + hostID]
        nodeDC <- which(v$nodeparents == 2 * p$obs - 1 + hostID)
        v$nodeparents[nodeDC] <- nodeUC
        # second, place it between the proposed upstream and downstream nodes
        nodeUP <- c(.ptr(v$nodeparents, hostID)[v$nodetimes[.ptr(v$nodeparents, hostID)] < tinf.prop], 0)[1]
        nodeDP <-intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID))
        v$nodeparents[nodeDP] <- 2 * p$obs - 1 + hostID
        v$nodeparents[2 * p$obs - 1 + hostID] <- nodeUP
        
        # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
        nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostioriorID)
        for (i in nodesinvolved) {
          if (any(.ptr(v$nodeparents, i) == 2 * p$obs - 1 + hostID)) {
            v$nodehosts[i] <- hostID
          } else {
            v$nodehosts[i] <- hostioriorID
          }
        }
        v$nodehosts[2 * p$obs - 1 + hostID] <- hostioriorID
        
        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate proposal ratio
        logproposalratio <- dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[hostID + 2 * p$obs - 1], 
                                   shape = tinf.prop.shape.mult * p$shape.sample, 
                                   scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) + 
          dgamma(.pbe0$v$nodetimes[hostiorID] - .pbe0$v$nodetimes[hostiorID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostiorID] - v$nodetimes[hostiorID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE)

        
        ### calculate likelihood
        .propose.pbe("trans")
        
        ### calculate acceptance probability
        logaccprob <- .pbe1$logLikgen + .pbe1$logLiksam + .pbe1$logLikcoal - .pbe0$logLikgen - .pbe0$logLiksam - .pbe0$logLikcoal + 
            logproposalratio
        
        ### accept or reject
        if (runif(1) < exp(logaccprob)) {
            .accept.pbe("trans")
        }
    }
    
    
}
