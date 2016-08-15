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
    p <- .pbe1$p
    v <- .pbe1$v
    
    ### propose the new infection time
    tinf.prop <- v$nodetimes[hostID] - 
      rgamma(1, shape = tinf.prop.shape.mult * p$shape.sample, scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample))
    .copy2pbe1("tinf.prop", le)
    
    ### going down the decision tree
    if (v$nodehosts[hostID + 2 * p$obs - 1] == 0) {
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
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        tinf.prop <- .pbe1$tinf.prop
        
        ### change to proposal state
            
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
            
        # phylotree in hostID
        v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <- 
          v$nodetimes[hostID + 2 * p$obs - 1] + 
          .samplecoaltimes(v$nodetimes[v$nodehosts == hostID & v$nodetypes != "c"] - v$nodetimes[hostID + 2 * p$obs - 1], 
                           p$wh.model, p$wh.slope)
        v$nodeparents[v$nodehosts == hostID] <- 
          .sampletopology(which(v$nodehosts == hostID), 
                          v$nodetimes[v$nodehosts == hostID], 
                          v$nodetypes[v$nodehosts == hostID], 
                          hostID + 2 * p$obs - 1, p$wh.model)
        
        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate proposal ratio
        logproposalratio <- dgamma(.pbe0$v$nodetimes[.pbe1$hostID] - .pbe0$v$nodetimes[.pbe1$hostID + 2 * p$obs - 1], 
                                   shape = tinf.prop.shape.mult * p$shape.sample, 
                                   scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1], 
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
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        tinf.prop <- .pbe1$tinf.prop

        ### change to proposal state
            
        # identify the current first infectee of hostID
        iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
        infectee.current.ID <- iees.nodeIDs[v$nodetimes[iees.nodeIDs] == min(v$nodetimes[iees.nodeIDs])] - 2 * p$obs + 1
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
            
        # propose infector for hostID
        dens.infectorproposal <- dgamma(tinf.prop - tail(v$nodetimes, p$obs), 
                                         shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) + 
                (tinf.prop - tail(v$nodetimes, p$obs) > 0)/.pbe0$h$dist[hostID, ]
        dens.infectorproposal[hostID] <- 0
        infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
            
        # bookkeeping: change the transmission tree topology
        v$nodehosts[infectee.current.ID + 2 * p$obs - 1] <- 0
        v$nodeparents[infectee.current.ID + 2 * p$obs - 1] <- 0
        v$nodehosts[hostID + 2 * p$obs - 1] <- infector.proposed.ID
            
        # prepare for phylotree proposals by moving the coalescent node from hostID to the new infector
        v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infector.proposed.ID
            
        # propose phylotrees for hostID and the new infector
        for (ID in c(hostID, infector.proposed.ID)) {
            v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[ID + 2 * p$obs - 1] + 
              .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2 * p$obs - 1], 
                               p$wh.model, p$wh.slope)
            v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], 
              v$nodetypes[v$nodehosts == ID], ID + 2 * p$obs - 1, p$wh.model)
        }
        
        ### update proposal environment
        .copy2pbe1("v", le)

        ### calculate proposal ratio
        logproposalratio <- log(sum(dens.infectorproposal)/(dens.infectorproposal[infector.proposed.ID])) + 
          dgamma(.pbe0$v$nodetimes[.pbe1$hostID] - .pbe0$v$nodetimes[.pbe1$hostID + 2 * p$obs - 1],
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1],
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
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        tinf.prop <- .pbe1$tinf.prop
      
        ### change to proposal state
            
        # identify the current first infectee of hostID
        iees.nodeIDs <- which(v$nodehosts == hostID & v$nodetypes == "t")
        infectee.current.ID <- iees.nodeIDs[v$nodetimes[iees.nodeIDs] == min(v$nodetimes[iees.nodeIDs])] - 2 * p$obs + 1
            
        # identify the other infectees of hostID and those of infectee.current.ID
        infecteenodes.primary <- iees.nodeIDs[iees.nodeIDs != (infectee.current.ID + 2 * p$obs - 1)]
        infecteenodes.secondary <- which(v$nodehosts == infectee.current.ID & v$nodetypes == "t")
            
        # remember current infector and infection time of hostID
        infector.current.ID <- v$nodehosts[hostID + 2 * p$obs - 1]
        inftime.current.hostID <- v$nodetimes[hostID + 2 * p$obs - 1]
            
        # bookkeeping: change infection times
        v$nodetimes[hostID + 2 * p$obs - 1] <- v$nodetimes[infectee.current.ID + 2 * p$obs - 1]
        v$nodetimes[infectee.current.ID + 2 * p$obs - 1] <- inftime.current.hostID
            
        # bookkeeping: change the transmission tree topology
        v$nodehosts[hostID + 2 * p$obs - 1] <- infectee.current.ID
        v$nodehosts[infectee.current.ID + 2 * p$obs - 1] <- infector.current.ID
        v$nodeparents[infectee.current.ID + 2 * p$obs - 1] <- v$nodeparents[hostID + 2 * p$obs - 1]
        if (exchange) {
            v$nodehosts[infecteenodes.primary] <- infectee.current.ID
            v$nodehosts[infecteenodes.secondary] <- hostID
        }
            
        # prepare for phylotree proposals by moving coalescent nodes between hostID and infectee.current.ID
        if (exchange) {
            coalnodes.current.ID <- v$nodehosts == hostID & v$nodetypes == "c"
            v$nodehosts[v$nodehosts == infectee.current.ID & v$nodetypes == "c"] <- hostID
            v$nodehosts[coalnodes.current.ID] <- infectee.current.ID
        } else {
            v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infectee.current.ID
        }
          
        # propose phylotrees for hostID and the new infector
        for (ID in c(hostID, infectee.current.ID)) {
            v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[ID + 2 * p$obs - 1] +
              .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2 * p$obs - 1], 
                               p$wh.model, p$wh.slope)
            v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], 
              v$nodetypes[v$nodehosts == ID], ID + 2 * p$obs - 1, p$wh.model)
        }

        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate proposal ratio
        if (infector.current.ID != 0) {
            ## .updatepathF: probability to sample after FIRST infectee
            logproposalratio <- pgamma(v$nodetimes[infectee.current.ID] - v$nodetimes[hostID + 2 * p$obs - 1],
                                       shape = tinf.prop.shape.mult * p$shape.sample, 
                                       scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log.p = TRUE) - 
              pgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1],
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
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        tinf.prop <- .pbe1$tinf.prop
      

        ### change to proposal state
        
        # identify the current index and current infector
        index.current.ID <- which(v$nodehosts == 0) - 2 * p$obs + 1
        infector.current.ID <- v$nodehosts[hostID + 2 * p$obs - 1]
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
        
        # bookkeeping: change the transmission tree topology
        v$nodehosts[hostID + 2 * p$obs - 1] <- 0
        v$nodeparents[hostID + 2 * p$obs - 1] <- 0
        v$nodehosts[index.current.ID + 2 * p$obs - 1] <- hostID
        
        # prepare for phylotree proposals by moving coalescent nodes from current infector to hostID
        v$nodehosts[v$nodehosts == infector.current.ID & v$nodetypes == "c"][1] <- hostID
        
        # propose phylotrees for hostID and the new infector
        for (ID in c(hostID, infector.current.ID)) {
          v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[ID + 2 * p$obs - 1] + 
            .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2 * p$obs - 1], 
                             p$wh.model, p$wh.slope)
          v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID],
                                                              v$nodetypes[v$nodehosts == ID], ID + 2 * p$obs - 1, p$wh.model)
        }

        ### update proposal environment
        .copy2pbe1("v", le)
        
        ### calculate proposal ratio the reverse proposal includes proposing an infector
        dens.infectorcurrent <- dgamma(.pbe0$v$nodetimes[hostID + 2 * p$obs - 1] - tail(.pbe0$v$nodetimes, p$obs), 
                                       shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) + 
          (.pbe0$v$nodetimes[hostID + 2 * p$obs - 1] - tail(.pbe0$v$nodetimes, p$obs) > 0)/h$dist[hostID, ]
        dens.infectorcurrent[hostID] <- 0
        
        logproposalratio <- log(dens.infectorcurrent[infector.current.ID]/(sum(dens.infectorcurrent))) + 
          dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[hostID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1], 
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
        h <- .pbe0$h
        p <- .pbe1$p
        v <- .pbe1$v
        hostID <- .pbe1$hostID
        tinf.prop <- .pbe1$tinf.prop
      

        ### change to proposal state
        
        # identify the current infector
        infector.current.ID <- v$nodehosts[hostID + 2 * p$obs - 1]
        
        # bookkeeping: change infection time
        v$nodetimes[hostID + 2 * p$obs - 1] <- tinf.prop
        
        # propose infector for hostID
        dens.infectorproposal <- dgamma(tinf.prop - tail(v$nodetimes, p$obs), 
                                        shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) +
          (tinf.prop - tail(v$nodetimes, p$obs) > 0)/h$dist[hostID, ]
        dens.infectorproposal[hostID] <- 0
        infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
        
        # bookkeeping: change the transmission tree topology
        v$nodehosts[hostID + 2 * p$obs - 1] <- infector.proposed.ID
        
        # prepare for phylotree proposals by moving the coalescent node from the current infector to the proposed infector
        v$nodehosts[v$nodehosts == infector.current.ID & v$nodetypes == "c"][1] <- infector.proposed.ID
        
        # propose phylotrees for hostID and the old and new infectors
        for (ID in unique(c(hostID, infector.current.ID, infector.proposed.ID))) {
          v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <- v$nodetimes[ID + 2 * p$obs - 1] + 
            .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2 * p$obs - 1], 
                             p$wh.model, p$wh.slope)
          v$nodeparents[v$nodehosts == ID] <- .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], 
                                                              v$nodetypes[v$nodehosts == ID], ID + 2 * p$obs - 1, p$wh.model)
        }

        ### update proposal environment
        .copy2pbe1("v", le)

        ### calculate proposal ratio the reverse proposal includes proposing an infector
        dens.infectorcurrent <- dgamma(.pbe0$v$nodetimes[hostID + 2 * p$obs - 1] - tail(.pbe0$v$nodetimes, p$obs), 
                                       shape = p$shape.gen, scale = p$mean.gen/p$shape.gen) +
          (.pbe0$v$nodetimes[hostID + 2 * p$obs - 1] - tail(.pbe0$v$nodetimes, p$obs) > 0)/h$dist[hostID, ]
        dens.infectorcurrent[hostID] <- 0
        
        logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal)/
                                  (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) + 
          dgamma(.pbe0$v$nodetimes[hostID] - .pbe0$v$nodetimes[hostID + 2 * p$obs - 1], 
                 shape = tinf.prop.shape.mult * p$shape.sample, 
                 scale = p$mean.sample/(tinf.prop.shape.mult * p$shape.sample), log = TRUE) - 
          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2 * p$obs - 1], shape = tinf.prop.shape.mult * p$shape.sample, 
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
