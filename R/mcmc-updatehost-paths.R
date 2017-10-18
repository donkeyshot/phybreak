### functions exclusively involved in updating the tree, i.e. proposing an infection time and infector, 
### and accepting or rejecting the proposal. The actual change in tree topology is done in the 'rewire_' functions.
### fixed parameters
tinf.prop.shape.mult <- 2/3  #shape for proposing infection time is sample.shape * tinf.prop.shape.mult

### updating both transmission and phylotree: proposing an infection time and following the decision tree
update_host <- function(hostID) {
    ### create an up-to-date proposal-environment with hostID as focal host
    prepare_pbe()
    copy2pbe1("hostID", environment())
  
    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    p <- pbe1$p
    v <- pbe1$v
    
    ### propose the new infection time
    tinf.prop <- v$nodetimes[hostID] - 
      rgamma(1, shape = tinf.prop.shape.mult * pbe1$p$sample.shape, scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape))
    copy2pbe1("tinf.prop", le)

    ### going down the decision tree
    if (v$infectors[hostID] == 0) {
      # Y (hostID is index case)
      if (tinf.prop < min(c(v$inftimes[v$infectors == hostID], Inf))) {
        # YY (... & tinf.prop before hostID's first transmission node)
        update_pathA()
      } else {
        # YN (... & tinf.prop after hostID's first transmission node)
        if (tinf.prop < sort(c(v$inftimes[v$infectors == hostID], Inf))[2]) {
          # YNY (... & tinf.prop before hostID's second transmission node)
          update_pathB()
        } else {
          # YNN (... & tinf.prop after hostID's second transmission node)
          if (runif(1) < 0.5) 
            update_pathC(TRUE) else update_pathC(FALSE)
        }
      }
    } else {
      # N (hostID is not index case)
      if (tinf.prop < v$inftimes[v$infectors == 0]) {
        # NY (... & tinf.prop before infection of index case)
        update_pathD()
      } else {
        # NN (... & tinf.prop after infection of index case)
        if (tinf.prop < min(c(v$inftimes[v$infectors == hostID], Inf))) {
          # NNY (... & tinf.prop before hostID's first transmission node)
          update_pathE()
        } else {
          # NNN (... & tinf.prop after hostID's first transmission node)
          if (runif(1) < 0.5) 
            update_pathF(TRUE) else update_pathF(FALSE)
        }
      }
    }
}


### updating transmission tree but keeping phylotree: proposing an infection time and following the decision tree 
update_host_keepphylo <- function(hostID) {
  ### create an up-to-date proposal-environment with hostID as focal host
  prepare_pbe()
  copy2pbe1("hostID", environment())
  
  ### making variables and parameters available within the function
  le <- environment()
  p <- pbe1$p
  v <- pbe1$v
  
  ### propose the new infection time
  tinf.prop <- v$nodetimes[hostID] - 
    rgamma(1, shape = tinf.prop.shape.mult * pbe1$p$sample.shape, scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape))
  copy2pbe1("tinf.prop", le)
  
  ### identify the focal host's infector
  hostiorID <- v$infectors[hostID]
  copy2pbe1("hostiorID", le)
  
  ### going down the decision tree
  if (hostiorID == 0) {
    # Y (hostID is index)
    if (tinf.prop < v$nodetimes[v$nodehosts == 0]) {
      # YY (... & tinf.prop before root of phylotree)
      update_pathG()
    } else {
      # YN (... & tinf.prop after root of phylotree): reject
    }
  } else {
    # N (hostID is not index)
    nodemrca <- .ptr(v$nodeparents, hostID)[.ptr(v$nodeparents, hostID) %in% .ptr(v$nodeparents, hostiorID)][1]
    timemrca <- v$nodetimes[nodemrca]
    copy2pbe1("nodemrca", le)
    copy2pbe1("timemrca", le)
    if (tinf.prop > timemrca) {
      # NY (... & tinf.prop after MRCA of hostID and hostiorID)
      update_pathH()
    } else {
      # NN (... & tinf.prop before MRCA of hostID and hostiorID)
      tinf2.prop <- v$nodetimes[hostiorID] - 
        rgamma(1, shape = tinf.prop.shape.mult * p$sample.shape, scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape))
      copy2pbe1("tinf2.prop", le)
      if (tinf2.prop > timemrca) {
        # NNY (... & tinf2.prop after MRCA of hostID and hostiorID)
        hostioriorID <- v$infectors[hostiorID]
        copy2pbe1("hostioriorID", le)
        if (hostioriorID == 0) {
          # NNYY (... & hostiorID is index)
          tinf.prop <- v$inftimes[hostiorID]
          copy2pbe1("tinf.prop", le)
          update_pathI()
        } else {
          # NNYN (... & hostiorID is not index)
          nodemrcaior <- .ptr(v$nodeparents, hostID)[.ptr(v$nodeparents, hostID) %in% .ptr(v$nodeparents, hostioriorID)][1]
          timemrcaior <- v$nodetimes[nodemrcaior]
          copy2pbe1("nodemrcaior", le)
          copy2pbe1("timemrcaior", le)
          if (tinf.prop > timemrcaior) {
            # NNYNY (... & tinf.prop after MRCA of hostID and hostioriorID)
            update_pathJ()
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


### updating one phylogenetic minitree while keeping the transmission tree
update_host_withinhost <- function(hostID) {
  ### create an up-to-date proposal-environment with hostID as focal host
  prepare_pbe()
  copy2pbe1("hostID", environment())
  
  ### change the tree
  rewire_pathK(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))
  
  ### calculate proposal ratio
  logproposalratio <- 0
  
  ### calculate likelihood
  propose_pbe("withinhost")
  
  ### calculate acceptance probability
  logaccprob <- pbe1$logLikseq - pbe0$logLikseq + logproposalratio
  
  ### accept or reject
  if (runif(1) < exp(logaccprob)) {
    accept_pbe("withinhost")
  }
}


{
  ### update if hostID is index and tinf.prop is before the first secondary case
  update_pathA <- function() {
    rewire_pathA(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))

    ### calculate proposal ratio
    logproposalratio <- dgamma(pbe0$v$nodetimes[pbe1$hostID] - pbe0$v$inftimes[pbe1$hostID],
                               shape = tinf.prop.shape.mult * pbe0$p$sample.shape,
                               scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape), log = TRUE) -
      dgamma(pbe1$v$nodetimes[pbe1$hostID] - pbe1$v$inftimes[pbe1$hostID],
             shape = tinf.prop.shape.mult * pbe1$p$sample.shape,
             scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape), log = TRUE)
    

    ### calculate likelihood
    propose_pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam - pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is index and tinf.prop is after the first secondary case, but before the second secondary
  update_pathB <- function() {
    ### propose infector for hostID
    dens.infectorproposal <- dgamma(pbe1$tinf.prop - pbe1$v$inftimes,
                                    shape = pbe1$p$gen.shape, scale = pbe1$p$gen.mean/pbe1$p$gen.shape) +
      (pbe1$tinf.prop - pbe1$v$inftimes > 0)/pbe1$h$dist[pbe1$hostID, ]
    dens.infectorproposal[pbe1$hostID] <- 0
    infector.proposed.ID <- sample(pbe1$p$obs, 1, prob = dens.infectorproposal)
    copy2pbe1("infector.proposed.ID", environment())

    ### change the trees
    rewire_pathB(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))

    ### calculate proposal ratio
    logproposalratio <- log(sum(dens.infectorproposal)/(dens.infectorproposal[infector.proposed.ID])) +
      dgamma(pbe0$v$nodetimes[pbe1$hostID] - pbe0$v$inftimes[pbe1$hostID],
             shape = tinf.prop.shape.mult * pbe0$p$sample.shape,
             scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape), log = TRUE) -
      dgamma(pbe1$v$nodetimes[pbe1$hostID] - pbe1$v$inftimes[pbe1$hostID],
             shape = tinf.prop.shape.mult * pbe1$p$sample.shape,
             scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape), log = TRUE)
    
    ### calculate likelihood
    propose_pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam - pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is index and tinf.prop is after the second secondary case 
  update_pathC <- function(exchange) {
    ### change the trees
    if(exchange) {
      rewire_pathCF2(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))
    } else {
      rewire_pathCF1(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))
    }
    
    ### calculate proposal ratio
    # second infectee not necessarily the same for reversal proposal, so identify these intervals for hostID and new index
    infectees.hostID <- which(pbe0$v$infectors == pbe1$hostID)
    infectee.first.ID <- infectees.hostID[pbe0$v$inftimes[infectees.hostID] == min(pbe0$v$inftimes[infectees.hostID])]
    sampleinterval.hostID <- pbe0$v$nodetimes[pbe1$hostID] - sort(pbe0$v$inftimes[infectees.hostID])[2]
    infectees.newindex <- which(pbe1$v$infectors == infectee.first.ID)
    sampleinterval.newindex <- pbe1$v$nodetimes[infectee.first.ID] - sort(c(pbe1$v$inftimes[infectees.newindex], Inf))[2]
    
    logproposalratio <- pgamma(sampleinterval.newindex, shape = tinf.prop.shape.mult * pbe0$p$sample.shape,
                               scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape), log.p = TRUE) -
      pgamma(sampleinterval.hostID, shape = tinf.prop.shape.mult * pbe1$p$sample.shape,
             scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape), log.p = TRUE)

    ### calculate likelihood
    propose_pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam - pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam +
      logproposalratio
    
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is before infection of the index case
  update_pathD <- function() {
    ### change the trees
    rewire_pathD(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))

    ### identify the current index and current infector
    index.current.ID <- which(pbe0$v$infectors == 0)
    infector.current.ID <- pbe0$v$infectors[pbe1$hostID]
    
    ### calculate proposal ratio 
    # the reverse proposal includes proposing an infector
    dens.infectorcurrent <- dgamma(pbe0$v$inftimes[pbe1$hostID] - pbe0$v$inftimes,
                                   shape = pbe0$p$gen.shape, scale = pbe0$p$gen.mean/pbe0$p$gen.shape) +
      (pbe0$v$inftimes[pbe1$hostID] - pbe0$v$inftimes > 0)/pbe0$h$dist[pbe1$hostID, ]
    dens.infectorcurrent[pbe1$hostID] <- 0
    
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID]/(sum(dens.infectorcurrent))) +
      dgamma(pbe0$v$nodetimes[pbe1$hostID] - pbe0$v$inftimes[pbe1$hostID],
             shape = tinf.prop.shape.mult * pbe0$p$sample.shape,
             scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape), log = TRUE) -
      dgamma(pbe1$v$nodetimes[pbe1$hostID] - pbe1$v$inftimes[pbe1$hostID],
             shape = tinf.prop.shape.mult * pbe1$p$sample.shape,
             scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape), log = TRUE)
    
    ### calculate likelihood
    propose_pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam - pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is after infection of the index case, but before the first secondary case
  update_pathE <- function() {
    ### identify the current infector and propose the new infector
    infector.current.ID <- pbe0$v$infectors[pbe1$hostID]
    dens.infectorproposal <- dgamma(pbe1$tinf.prop - pbe1$v$inftimes,
                                    shape = pbe1$p$gen.shape, scale = pbe1$p$gen.mean/pbe1$p$gen.shape) +
      (pbe1$tinf.prop - pbe1$v$inftimes > 0)/pbe1$h$dist[pbe1$hostID, ]
    dens.infectorproposal[pbe1$hostID] <- 0
    infector.proposed.ID <- sample(pbe1$p$obs, 1, prob = dens.infectorproposal)
    copy2pbe1("infector.proposed.ID", environment())
    
    ### change the trees    
    rewire_pathE(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))
 
    ### calculate proposal ratio 
    # the reverse proposal includes proposing an infector
    dens.infectorcurrent <- dgamma(pbe0$v$inftimes[pbe1$hostID] - pbe0$v$inftimes,
                                   shape = pbe0$p$gen.shape, scale = pbe0$p$gen.mean/pbe0$p$gen.shape) +
      (pbe0$v$inftimes[pbe1$hostID] - pbe0$v$inftimes > 0)/pbe0$h$dist[pbe1$hostID, ]
    dens.infectorcurrent[pbe1$hostID] <- 0
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal)/
                              (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) +
      dgamma(pbe0$v$nodetimes[pbe1$hostID] - pbe0$v$inftimes[pbe1$hostID],
             shape = tinf.prop.shape.mult * pbe0$p$sample.shape,
             scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape), log = TRUE) -
      dgamma(pbe1$v$nodetimes[pbe1$hostID] - pbe1$v$inftimes[pbe1$hostID], 
             shape = tinf.prop.shape.mult * pbe1$p$sample.shape,
             scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape), log = TRUE)
    
    
    ### calculate likelihood
    propose_pbe("phylotrans")

    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam - pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("phylotrans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is after the first secondary case
  update_pathF <- function(exchange) {
      ### change the trees
      if(exchange) {
      rewire_pathCF2(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))
    } else {
      rewire_pathCF1(pbe0$p$wh.model %in% c(4, 5, "exponential", "constant"))
    }
    
    ### calculate proposal ratio
    infectees.hostID <- which(pbe0$v$infectors == pbe1$hostID)
    infectee.first.ID <- infectees.hostID[pbe0$v$inftimes[infectees.hostID] == min(pbe0$v$inftimes[infectees.hostID])]
    logproposalratio <- pgamma(pbe0$v$nodetimes[infectee.first.ID] - pbe0$v$inftimes[pbe1$hostID],
                               shape = tinf.prop.shape.mult * pbe0$p$sample.shape,
                               scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape), log.p = TRUE) -
      pgamma(pbe1$v$nodetimes[pbe1$hostID] - pbe1$v$inftimes[infectee.first.ID],
             shape = tinf.prop.shape.mult * pbe1$p$sample.shape,
             scale = pbe1$p$sample.mean/(tinf.prop.shape.mult * pbe1$p$sample.shape), log.p = TRUE)
    
    ### calculate likelihood
    propose_pbe("phylotrans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam - pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam +
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("phylotrans")
    }
  }
  
  
}



{
  
  ### update if hostID is index and tinf.prop is before the first coalescence node
  update_pathG <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop
    
    
    ### change to proposal state
    
    # bookkeeping: change infection time
    v$nodetimes[hostID + 2 * d$nsamples - 1] <- tinf.prop
    v$inftimes[hostID] <- tinf.prop
    
    
    ### calculate proposal ratio
    logproposalratio <- dgamma(pbe0$v$nodetimes[hostID] - pbe0$v$inftimes[hostID], 
                               shape = tinf.prop.shape.mult * p$sample.shape, 
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) - 
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID], 
             shape = tinf.prop.shape.mult * p$sample.shape, 
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    
    ### update proposal environment
    copy2pbe1("v", le)
    
    ### calculate likelihood
    propose_pbe("trans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen + pbe1$logLiksam + pbe1$logLikcoal - pbe0$logLikgen - pbe0$logLiksam - pbe0$logLikcoal + 
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("trans")
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is after the MRCA of hostID and infector 
  update_pathH <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    hostID <- pbe1$hostID
    hostiorID <- pbe1$hostiorID
    tinf.prop <- pbe1$tinf.prop
    
    ### change to proposal state
    
    # bookkeeping: change infection time
    v$nodetimes[hostID + 2 * d$nsamples - 1] <- tinf.prop
    v$inftimes[hostID] <- tinf.prop
    
    # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
    nodeUC <- v$nodeparents[2 * d$nsamples - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2 * d$nsamples - 1 + hostID)
    v$nodeparents[nodeDC] <- nodeUC
    # second, place it between the proposed upstream and downstream nodes
    nodeUP <- .ptr(v$nodeparents, hostID)[v$nodetimes[.ptr(v$nodeparents, hostID)] < tinf.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID))
    v$nodeparents[nodeDP] <- 2 * d$nsamples - 1 + hostID
    v$nodeparents[2 * d$nsamples - 1 + hostID] <- nodeUP
    
    # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostiorID)
    for (i in nodesinvolved) {
      if (any(.ptr(v$nodeparents, i) == 2 * d$nsamples - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostiorID
      }
    }
    v$nodehosts[2 * p$obs - 1 + hostID] <- hostiorID
    v$infectors <- tail(v$nodehosts, p$obs)
    
    ### update proposal environment
    copy2pbe1("v", le)
    
    ### calculate proposal ratio
    logproposalratio <- dgamma(pbe0$v$nodetimes[hostID] - pbe0$v$inftimes[hostID],
                               shape = tinf.prop.shape.mult * p$sample.shape, 
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) - 
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID], 
             shape = tinf.prop.shape.mult * p$sample.shape, 
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    
    
    ### calculate likelihood
    propose_pbe("trans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen + pbe1$logLiksam + pbe1$logLikcoal - pbe0$logLikgen - pbe0$logLiksam - pbe0$logLikcoal + 
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("trans")
    }
  }
  
  
  ### update if hostID is not index, tinf.prop is before the MRCA of hostID and infector, and infector is index 
  update_pathI <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    hostID <- pbe1$hostID
    hostiorID <- pbe1$hostiorID
    hostioriorID <- pbe1$hostioriorID
    tinf.prop <- pbe1$tinf.prop
    tinf2.prop <- pbe1$tinf2.prop
    timemrca <- pbe1$timemrca
    
    ### change to proposal state
    
    ## first the infector: move transmission node downstream
    
    # bookkeeping: change infection time
    v$nodetimes[hostiorID + 2 * p$obs - 1] <- tinf2.prop
    v$inftimes[hostiorID] <- tinf2.prop
    
    # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
    nodeUC <- v$nodeparents[2 * d$nsamples - 1 + hostiorID]
    nodeDC <- which(v$nodeparents == 2 * d$nsamples - 1 + hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    # second, place it between the proposed upstream and downstream nodes
    nodeUP <- .ptr(v$nodeparents, hostiorID)[v$nodetimes[.ptr(v$nodeparents, hostiorID)] < tinf2.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostiorID))
    v$nodeparents[nodeDP] <- 2 * d$nsamples - 1 + hostiorID
    v$nodeparents[2 * d$nsamples - 1 + hostiorID] <- nodeUP
    
    # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
    nodesinvolved <- which(v$nodehosts == hostiorID | v$nodehosts == hostioriorID)
    for (i in nodesinvolved) {
      if (any(.ptr(v$nodeparents, i) == 2 * d$nsamples - 1 + hostiorID)) {
        v$nodehosts[i] <- hostiorID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2 * d$nsamples - 1 + hostiorID] <- hostioriorID
    
    ## then hostID itself: move transmission node upstream
    
    # bookkeeping: change infection time
    v$nodetimes[hostID + 2 * d$nsamples - 1] <- tinf.prop
    v$inftimes[hostID] <- tinf.prop
    
    # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
    nodeUC <- v$nodeparents[2 * d$nsamples - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2 * d$nsamples - 1 + hostID)
    v$nodeparents[nodeDC] <- nodeUC
    # second, place it between the proposed upstream and downstream nodes
    nodeUP <- c(.ptr(v$nodeparents, hostID)[v$nodetimes[.ptr(v$nodeparents, hostID)] < tinf.prop], 0)[1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID))
    v$nodeparents[nodeDP] <- 2 * d$nsamples - 1 + hostID
    v$nodeparents[2 * d$nsamples - 1 + hostID] <- nodeUP
    
    # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostioriorID)
    for (i in nodesinvolved) {
      if (any(.ptr(v$nodeparents, i) == 2 * d$nsamples - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2 * d$nsamples - 1 + hostID] <- hostioriorID
    v$infectors <- tail(v$nodehosts, p$obs)
    
    ### update proposal environment
    copy2pbe1("v", le)
    
    ### calculate proposal ratio
    logproposalratio <- dgamma(pbe0$v$nodetimes[hostID] - pbe0$v$inftimes[hostID],
                               shape = tinf.prop.shape.mult * p$sample.shape,
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) - 
      log(1 - pgamma(v$nodetimes[hostID] - timemrca,
                     shape = tinf.prop.shape.mult * p$sample.shape, 
                     scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape))) + 
      log(1 - pgamma(pbe0$v$nodetimes[hostiorID] - timemrca, 
                     shape = tinf.prop.shape.mult * p$sample.shape, 
                     scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape))) - 
      dgamma(v$nodetimes[hostiorID] - v$inftimes[hostiorID], 
             shape = tinf.prop.shape.mult * p$sample.shape, 
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    
    
    
    ### calculate likelihood
    propose_pbe("trans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen + pbe1$logLiksam + pbe1$logLikcoal - pbe0$logLikgen - pbe0$logLiksam - pbe0$logLikcoal + 
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("trans")
    }
  }
  
  
  ### update if hostID is not index, tinf.prop is before the MRCA of hostID and infector, infector is not index, and tinf.prop is
  ### after the MRCA of hostID and infector's infector 
  update_pathJ <- function() {
    ### making variables and parameters available within the function
    le <- environment()
    d <- pbe0$d
    h <- pbe0$h
    p <- pbe1$p
    v <- pbe1$v
    hostID <- pbe1$hostID
    hostiorID <- pbe1$hostiorID
    hostioriorID <- pbe1$hostioriorID
    tinf.prop <- pbe1$tinf.prop
    tinf2.prop <- pbe1$tinf2.prop
    
    ### change to proposal state
    
    ## first the infector: move transmission node downstream
    
    # bookkeeping: change infection time
    v$nodetimes[hostiorID + 2 * d$nsamples - 1] <- tinf2.prop
    v$inftimes[hostiorID] <- tinf2.prop
    
    # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
    nodeUC <- v$nodeparents[2 * d$nsamples - 1 + hostiorID]
    nodeDC <- which(v$nodeparents == 2 * d$nsamples - 1 + hostiorID)
    v$nodeparents[nodeDC] <- nodeUC
    # second, place it between the proposed upstream and downstream nodes
    nodeUP <- .ptr(v$nodeparents, hostiorID)[v$nodetimes[.ptr(v$nodeparents, hostiorID)] < tinf2.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostiorID))
    v$nodeparents[nodeDP] <- 2 * d$nsamples - 1 + hostiorID
    v$nodeparents[2 * d$nsamples - 1 + hostiorID] <- nodeUP
    
    # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
    nodesinvolved <- which(v$nodehosts == hostiorID | v$nodehosts == hostioriorID)
    for (i in nodesinvolved) {
      if (any(.ptr(v$nodeparents, i) == 2 * d$nsamples - 1 + hostiorID)) {
        v$nodehosts[i] <- hostiorID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2 * d$nsamples - 1 + hostiorID] <- hostioriorID
    
    ## then hostID itself: move transmission node upstream
    
    # bookkeeping: change infection time
    v$nodetimes[hostID + 2 * d$nsamples - 1] <- tinf.prop
    v$inftimes[hostID] <- tinf.prop
    
    # bookkeeping: move the transmission node in the phylotree first, remove it by connecting the upstream and downstream nodes
    nodeUC <- v$nodeparents[2 * d$nsamples - 1 + hostID]
    nodeDC <- which(v$nodeparents == 2 * d$nsamples - 1 + hostID)
    v$nodeparents[nodeDC] <- nodeUC
    # second, place it between the proposed upstream and downstream nodes
    nodeUP <- c(.ptr(v$nodeparents, hostID)[v$nodetimes[.ptr(v$nodeparents, hostID)] < tinf.prop], 0)[1]
    nodeDP <-intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID))
    v$nodeparents[nodeDP] <- 2 * d$nsamples - 1 + hostID
    v$nodeparents[2 * d$nsamples - 1 + hostID] <- nodeUP
    
    # bookkeeping: give all infectees of hostID and hostiorID their correct infector according to the proposed transmission node
    nodesinvolved <- which(v$nodehosts == hostID | v$nodehosts == hostioriorID)
    for (i in nodesinvolved) {
      if (any(.ptr(v$nodeparents, i) == 2 * d$nsamples - 1 + hostID)) {
        v$nodehosts[i] <- hostID
      } else {
        v$nodehosts[i] <- hostioriorID
      }
    }
    v$nodehosts[2 * d$nsamples - 1 + hostID] <- hostioriorID
    v$infectors <- tail(v$nodehosts, p$obs)
    
    ### update proposal environment
    copy2pbe1("v", le)
    
    ### calculate proposal ratio
    logproposalratio <- dgamma(pbe0$v$nodetimes[hostID] - pbe0$v$inftimes[hostID], 
                               shape = tinf.prop.shape.mult * p$sample.shape, 
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) - 
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID], 
             shape = tinf.prop.shape.mult * p$sample.shape, 
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) + 
      dgamma(pbe0$v$nodetimes[hostiorID] - pbe0$v$inftimes[hostiorID], 
             shape = tinf.prop.shape.mult * p$sample.shape, 
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) - 
      dgamma(v$nodetimes[hostiorID] - v$inftimes[hostiorID], 
             shape = tinf.prop.shape.mult * p$sample.shape, 
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    
    
    ### calculate likelihood
    propose_pbe("trans")
    
    ### calculate acceptance probability
    logaccprob <- pbe1$logLikgen + pbe1$logLiksam + pbe1$logLikcoal - pbe0$logLikgen - pbe0$logLiksam - pbe0$logLikcoal + 
      logproposalratio
    
    ### accept or reject
    if (runif(1) < exp(logaccprob)) {
      accept_pbe("trans")
    }
  }
  
  
}
