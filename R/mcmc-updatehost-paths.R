### functions exclusively involved in updating the tree, i.e. proposing an infection time and infector, 
### and accepting or rejecting the proposal. The actual change in tree topology is done in the 'rewire_' functions.
### fixed parameters
tinf.prop.shape.mult <- 2/3  #shape for proposing infection time is sample.shape * tinf.prop.shape.mult

### fork to the requested update protocol
update_host <- function(hostID, which_protocol, history) {
  ### use protocol
  if(which_protocol == "keepphylo") {
    update_host_keepphylo(hostID)
  } else if(which_protocol == "withinhost") {
    update_host_withinhost(hostID)
  } else {
    if(history)
      update_host_history(hostID, which_protocol)
    else 
      update_host_phylotrans(hostID, which_protocol)
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
  if(pbe0$p$wh.bottleneck == "wide") {
    rewire_pathK_wide_classic()
  } else {
    rewire_pathK_complete_classic()
  }

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


### updating both transmission and phylotree: proposing an infection time and following the decision tree
update_host_phylotrans <- function(hostID, which_protocol) { 
  ### copy focal host "hostID" to proposal-environment
  copy2pbe1("hostID", environment())
  
  ### making variables and parameters available within the function
  le <- environment()
  d <- pbe0$d
  p <- pbe0$p
  v <- pbe0$v
  
  ### propose the new infection time
  tinf.prop <- v$nodetimes[hostID] -
    rgamma(1, shape = tinf.prop.shape.mult * pbe0$p$sample.shape, scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape))
  # tinf.prop <- v$inftimes[hostID] + rnorm(1, 0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape))
  # tinf.prop <- min(tinf.prop, 2 * v$nodetimes[hostID] - tinf.prop)
  copy2pbe1("tinf.prop", le)
  
  ### going down the decision tree
  if (v$infectors[hostID] == 0) {
    # Y (hostID is index case)
    if (sum(v$infectors == hostID) == 0) {
      # YY (... & no transmission from hostID)
      update_pathA(which_protocol)
    } else{
      # YN (... & hostID not only case in subtree)
      if (tinf.prop < min(c(v$inftimes[v$infectors == hostID], Inf))) {
        # YNY (... & tinf.prop before hostID's first transmission node)
        update_pathA(which_protocol)
      } else {
        # YNN (... & tinf.prop after hostID's first transmission node)
        if (tinf.prop < sort(c(v$inftimes[v$infectors == hostID], Inf))[2]) {
          # YNNY (... & tinf.prop before hostID's second transmission node)
          update_pathB(which_protocol)
        } else {
          # YNNN (... & tinf.prop after hostID's second transmission node)
          update_pathC(which_protocol)
        }
      }
    }
  } else {
    # N (hostID is not index case)
    if (tinf.prop < v$inftimes[v$tree[hostID]]) {
      # NY (... & tinf.prop before infection of hostID's index case)
      update_pathD(which_protocol)
      # return()
    } else {
      # NN (... & tinf.prop after infection of hostID's index case)
      if (tinf.prop < min(c(v$inftimes[v$infectors == hostID], Inf))) {
        # NNY (... & tinf.prop before hostID's first transmission node)
        update_pathE(which_protocol)
      } else {
        # NNN (... & tinf.prop after hostID's first transmission node)
        update_pathF(which_protocol) 
      }
    }
  }
}


### updating index cases: proposing an infection time and add or remove index case
update_host_history <- function(hostID, which_protocol) {
  ### create an up-to-date proposal-environment with hostID as focal host
  copy2pbe1("hostID", environment())
  
  ### making variables and parameters available within the function
  le <- environment()
  d <- pbe0$d
  p <- pbe0$p
  v <- pbe0$v
  
  ### propose the new infection time
  tinf.prop <- v$nodetimes[hostID] -
    rgamma(1, shape = tinf.prop.shape.mult * pbe0$p$sample.shape, scale = pbe0$p$sample.mean/(tinf.prop.shape.mult * pbe0$p$sample.shape))
  copy2pbe1("tinf.prop", le)
  
  ### going down the decision tree
  if (hostID == 0) {
    # Y (hostID is history)
    update_host_withinhost(hostID)
  } else {
    # N (hostID is not history)
    if (tinf.prop < min(c(v$inftimes[v$infectors == hostID], Inf))) {
      # NY (... & tinf.prop before hostID's first transmission node)
      update_pathL(which_protocol)
    }
  }
}

{
  ### update if hostID is index and tinf.prop is before the first secondary case
  update_pathA <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    
    tinf.prop <- pbe1$tinf.prop
    
    ### calculate proposal ratio
    # logproposalratio <- 0
    logproposalratio <- dgamma(v$nodetimes[hostID] - v$inftimes[hostID],
                               shape = tinf.prop.shape.mult * p$sample.shape,
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) -
      dgamma(v$nodetimes[hostID] - tinf.prop,
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(which_protocol == "classic") {
      if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathA_complete_classic"
      } else {
        rewire_function <- "rewire_pathA_wide_classic"
      }
    } else if(p$wh.bottleneck == "complete") {
      rewire_function <- "rewire_pathA_complete_edgewise"
    } else {
      rewire_function <- "rewire_pathA_wide_edgewise"
    }
    update_move(rewire_function, which_protocol)
  }
  
  
  ### update if hostID is index and tinf.prop is after the first secondary case, but before the second secondary
  update_pathB <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop
    
    ### propose infector for hostID
    dens.infectorproposal <- infect_distribution(tinf.prop, 
                                                 v$inftimes, p,
                                                 nodetimes = v$nodetimes[v$nodetypes=="s"],
                                                 cultimes = v$cultimes) +
      (tinf.prop - v$inftimes > 0)/pbe0$h$dist[hostID, ]
    dens.infectorproposal[which(v$tree != v$tree[hostID])] <- 0
    dens.infectorproposal[hostID] <- 0
    infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
    copy2pbe1("infector.proposed.ID", environment())

    ### calculate proposal ratio
    # logproposalratio <- log(sum(dens.infectorproposal)/(dens.infectorproposal[infector.proposed.ID])) 
    logproposalratio <- log(sum(dens.infectorproposal)/(dens.infectorproposal[infector.proposed.ID])) +
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID],
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) -
      dgamma(v$nodetimes[hostID] - tinf.prop,
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(which_protocol == "classic") {
      if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathB_complete_classic"
      } else {
        rewire_function <- "rewire_pathB_wide_classic"
      }
    } else if(p$wh.bottleneck == "complete") {
      rewire_function <- "rewire_pathB_complete_edgewise"
    } else {
      rewire_function <- "rewire_pathB_wide_edgewise"
    }
    update_move(rewire_function, which_protocol)
  }
  
  
  ### update if hostID is index and tinf.prop is after the second secondary case 
  update_pathC <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop
    
    ### calculate proposal ratio
    # second infectee not necessarily the same for reversal proposal, so identify these intervals for hostID and new index
    infectees.hostID <- which(pbe0$v$infectors == hostID)
    sampleinterval.hostID <- pbe0$v$nodetimes[hostID] - sort(pbe0$v$inftimes[infectees.hostID])[2]
    newindexID <- infectees.hostID[pbe0$v$inftimes[infectees.hostID] == min(pbe0$v$inftimes[infectees.hostID])]
    sampleinterval.newindex <- if(which_protocol == "classic" && (exchange <- runif(1) < 0.5)) {
      sampleinterval.hostID
    } else {
      infectees.newindex <- which(pbe0$v$infectors == newindexID)
      pbe0$v$nodetimes[newindexID] - sort(c(pbe0$v$inftimes[infectees.newindex], Inf))[1]
    }

    # logproposalratio <- pnorm(v$inftimes[newindexID] - v$nodetimes[newindexID] + sampleinterval.newindex, 
    #                           0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape)) -
    #   pnorm(v$inftimes[newindexID] - v$nodetimes[newindexID] - sampleinterval.newindex, 
    #         0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape)) -
    #   pnorm(v$inftimes[hostID] - v$nodetimes[hostID] + sampleinterval.hostID, 
    #         0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape)) +
    #   pnorm(v$inftimes[hostID] - v$nodetimes[hostID] - sampleinterval.hostID, 
    #         0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape))
    logproposalratio <- pgamma(sampleinterval.newindex, shape = tinf.prop.shape.mult * p$sample.shape,
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log.p = TRUE) -
      pgamma(sampleinterval.hostID, shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log.p = TRUE)
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(logproposalratio > -Inf) {
      ### propose minitrees and accept or reject
      if(which_protocol == "classic") {
        if(p$wh.bottleneck == "complete") {
          rewire_function <- c("rewire_pathCF1_complete_classic", "rewire_pathCF2_complete_classic")[1 + exchange]
        } else {
          rewire_function <- c("rewire_pathCF1_wide_classic", "rewire_pathCF2_wide_classic")[1 + exchange]
        }
      } else if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathCF_complete_edgewise"
      } else {
        rewire_function <- "rewire_pathCF_wide_edgewise"
      }
      update_move(rewire_function, which_protocol)
    }
  }
  
  
  ### update if hostID is not index and tinf.prop is before infection of the index case
  update_pathD <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop
    
    ### calculate proposal ratio 
    # the reverse proposal includes proposing an infector, 
    # so first identify the current infector
    infector.current.ID <- v$infectors[hostID]
    dens.infectorcurrent <- infect_distribution(v$inftimes[hostID], 
                                                v$inftimes, p,
                                                nodetimes = v$nodetimes[v$nodetypes=="s"],
                                                cultimes = v$cultimes) +
      (v$inftimes[hostID] - v$inftimes > 0)/pbe0$h$dist[hostID, ]
    dens.infectorcurrent[which(v$tree != v$tree[hostID])] <- 0
    dens.infectorcurrent[hostID] <- 0
    
    # logproposalratio <- log(dens.infectorcurrent[infector.current.ID]/(sum(dens.infectorcurrent))) 
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID]/(sum(dens.infectorcurrent))) +
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID],
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) -
      dgamma(v$nodetimes[hostID] - tinf.prop,
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)    
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(which_protocol == "classic") {
      if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathD_complete_classic"
      } else {
        rewire_function <- "rewire_pathD_wide_classic"
      }
    } else if(p$wh.bottleneck == "complete") {
      rewire_function <- "rewire_pathD_complete_edgewise"
    } else {
      rewire_function <- "rewire_pathD_wide_edgewise"
    }
    update_move(rewire_function, which_protocol)
  }
  
  
  ### update if hostID is not index and tinf.prop is after infection of the index case, but before the first secondary case
  update_pathE <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop

    ### identify the current infector and propose the new infector
    infector.current.ID <- v$infectors[hostID]
    dens.infectorproposal <- infect_distribution(tinf.prop, 
                                                 v$inftimes, p,
                                                 nodetimes = v$nodetimes[v$nodetypes=="s"],
                                                 cultimes = v$cultimes) +
      (tinf.prop - v$inftimes > 0)/pbe0$h$dist[hostID, ]
    
    if(length(v$cultimes) > 0)
      dens.infectorproposal <- dens.infectorproposal * (v$cultimes + 5 - tinf.prop > 0)
  
    dens.infectorproposal[which(v$tree != v$tree[hostID])] <- 0
    dens.infectorproposal[hostID] <- 0
  
    infector.proposed.ID <- sample(p$obs, 1, prob = dens.infectorproposal)
    
    ### calculate proposal ratio 
    # the reverse proposal includes proposing an infector
    dens.infectorcurrent <- infect_distribution(v$inftimes[hostID], 
                                                v$inftimes, p,
                                                nodetimes = v$nodetimes[v$nodetypes=="s"],
                                                cultimes = v$cultimes) +
      (v$inftimes[hostID] - v$inftimes > 0)/pbe0$h$dist[hostID, ]

    dens.infectorcurrent[which(v$tree != v$tree[hostID])] <- 0
    dens.infectorcurrent[hostID] <- 0
    
    # logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal)/
    #                           (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) 
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal)/
                              (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) +
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID],
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) -
      dgamma(v$nodetimes[hostID] - tinf.prop,
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    
    copy2pbe1("infector.proposed.ID", environment())
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(which_protocol == "classic") {
      if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathE_complete_classic"
      } else {
        rewire_function <- "rewire_pathE_wide_classic"
      }
    } else if(p$wh.bottleneck == "complete") {
      rewire_function <- "rewire_pathE_complete_edgewise"
    } else {
      rewire_function <- "rewire_pathE_wide_edgewise"
    }
    update_move(rewire_function, which_protocol)
  }
  
  
  ### update if hostID is not index and tinf.prop is after the first secondary case
  update_pathF <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop
    
    
    ### calculate proposal ratio
    infectees.hostID <- which(v$infectors == hostID)
    infectee.first.ID <- infectees.hostID[v$inftimes[infectees.hostID] == min(v$inftimes[infectees.hostID])]
    # logproposalratio <-  -  pnorm(v$inftimes[hostID] - v$inftimes[infectee.first.ID] - 2 * v$nodetimes[infectee.first.ID], 
    #         0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape)) +
    #   pnorm(v$inftimes[hostID] - v$inftimes[infectee.first.ID] - 2 * v$nodetimes[hostID], 
    #         0, 0.5 * pbe0$h$mS.av / sqrt(p$sample.shape))
    logproposalratio <- pgamma(v$nodetimes[infectee.first.ID] - v$inftimes[infectee.first.ID],
                               shape = tinf.prop.shape.mult * p$sample.shape,
                               scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log.p = TRUE) -
      pgamma(v$nodetimes[hostID] - v$inftimes[infectee.first.ID],
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log.p = TRUE)
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(logproposalratio > -Inf) {
      ### propose minitrees and accept or reject
      if(which_protocol == "classic") {
        exchange <- runif(1) < 0.5
        if(p$wh.bottleneck == "complete") {
          rewire_function <- c("rewire_pathCF1_complete_classic", "rewire_pathCF2_complete_classic")[1 + exchange]
        } else {
          rewire_function <- c("rewire_pathCF1_wide_classic", "rewire_pathCF2_wide_classic")[1 + exchange]
        }
      } else if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathCF_complete_edgewise"
      } else {
        rewire_function <- "rewire_pathCF_wide_edgewise"
      }
      update_move(rewire_function, which_protocol)
    }
  }
  
  
  ### update if proposal is add/remove index case and tinf.prop if before first secondary case
  update_pathL <- function(which_protocol) {
    ### Make input locally available
    p <- pbe0$p
    v <- pbe0$v
    hostID <- pbe1$hostID
    tinf.prop <- pbe1$tinf.prop
    
    ### identify the current infector and propose the new infector
    infector.current.ID <- v$infectors[hostID]
    if (infector.current.ID == 0) infector.current.ID <- p$obs+1
    
    dens.infectorproposal <- c(infect_distribution(tinf.prop, 
                                                   v$inftimes, p,
                                                   nodetimes = v$nodetimes[v$nodetypes=="s"],
                                                   cultimes = v$cultimes) +
                                 (tinf.prop - v$inftimes > 0)/pbe0$h$dist[hostID, ], 1)
    
    if(length(v$cultimes) > 0)
      dens.infectorproposal <- dens.infectorproposal * (v$cultimes + 5 - tinf.prop > 0)
 
    dens.infectorproposal[hostID] <- 0
    
    infector.proposed.ID <- sample(p$obs+1, 1, prob = dens.infectorproposal)
    
    ### calculate proposal ratio 
    # the reverse proposal includes proposing an infector
    dens.infectorcurrent <- c(infect_distribution(v$inftimes[hostID], 
                                                  v$inftimes, p,
                                                  nodetimes = v$nodetimes[v$nodetypes=="s"],
                                                  cultimes = v$cultimes) +
                                (v$inftimes[hostID] - v$inftimes > 0)/pbe0$h$dist[hostID, ], 1)
    dens.infectorcurrent[hostID] <- 0
    
    logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal)/
                              (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) +
      dgamma(v$nodetimes[hostID] - v$inftimes[hostID],
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE) -
      dgamma(v$nodetimes[hostID] - tinf.prop,
             shape = tinf.prop.shape.mult * p$sample.shape,
             scale = p$sample.mean/(tinf.prop.shape.mult * p$sample.shape), log = TRUE)
    
    if (infector.proposed.ID == p$obs+1) infector.proposed.ID <- 0
    
    copy2pbe1("infector.proposed.ID", environment())
    copy2pbe1("logproposalratio", environment())
    
    ### propose minitrees and accept or reject
    if(which_protocol == "classic") {
      if(p$wh.bottleneck == "complete") {
        rewire_function <- "rewire_pathL_complete_classic"
      } else {
        rewire_function <- "rewire_pathL_wide_classic"
      }
    } else if(p$wh.bottleneck == "complete") {
      rewire_function <- "rewire_pathL_complete_edgewise"
    } else {
      rewire_function <- "rewire_pathL_wide_edgewise"
    }
    update_move(rewire_function, which_protocol)
  }
  
  ### update of phylogenetic tree
  update_move <- function(rewirefunction, which_protocol) {
    prepare_pbe()
    do.call(rewirefunction, args = list())
    
    if(pbe1$logLiktoporatio > -Inf) {
      propose_pbe("phylotrans")
      logacceptanceprob <- pbe0$heat * 
        (pbe1$logLikseq + pbe1$logLikgen + pbe1$logLiksam + pbe1$logLikdist + pbe1$logLiktoporatio -
        pbe0$logLikseq - pbe0$logLikgen - pbe0$logLiksam - pbe0$logLikdist) + pbe1$logproposalratio
      
      if (runif(1) < exp(logacceptanceprob)) {
        accept_pbe("phylotrans")
      }
    }
    
    if(which_protocol == "edgewise") {
      update_move_sampleedges()
    }
  }
  
  
  ### update after first step in two-step update protocol
  update_move_sampleedges <- function() {
    # sample edges in hostID
    sampleedges <- which(pbe0$v$nodehosts == pbe1$hostID & pbe0$v$nodetypes %in% c("s", "x"))
    sampleedges <- if(length(sampleedges) == 1) sampleedges else sample(sampleedges)
    
    # deconnect, reconnect, accept/reject each tip one by one
    for(edge in sampleedges) {
      if(pbe0$v$nodetypes[pbe0$v$nodeparents[edge]] == "c" || 
         (pbe0$p$wh.bottleneck == "wide" &&
          any(pbe0$v$nodetypes %in% c("x", "t") & pbe0$v$nodehosts == pbe1$hostID))) {
        prepare_pbe()
        if (!is.null(pbe1$infector.proposed.ID)) rm('infector.proposed.ID', envir = pbe1)
        coalnode <- take_cnode(edge)
        if(pbe1$v$nodehosts[2 * pbe1$d$nsamples - 1 + pbe1$hostID] == -1) {
          node_torelink <- which(pbe1$v$nodehosts == pbe1$hostID)
          node_torelink <- node_torelink[pbe1$v$nodetypes[pbe1$v$nodeparents[node_torelink]] == "b"][1]
          bnode_toremove <- pbe1$v$nodeparents[node_torelink]
          pbe1$v$nodeparents[c(node_torelink, 2 * pbe1$d$nsamples - 1 + pbe1$hostID, bnode_toremove)] <-
            c(2 * pbe1$d$nsamples - 1 + pbe1$hostID, pbe1$v$nodeparents[bnode_toremove], -1L)
          pbe1$v$nodehosts[c(2 * pbe1$d$nsamples - 1 + pbe1$hostID, bnode_toremove)] <- 
            c(pbe1$v$nodehosts[bnode_toremove], -1L)
          pbe1$v$nodetypes[c(2 * pbe1$d$nsamples - 1 + pbe1$hostID, bnode_toremove)] <- c("t", "0")
        }
        pbe1$v$nodehosts[edge] <- pbe1$hostID
        if(pbe0$p$wh.bottleneck == "wide") {
          rewire_pullnodes_wide(pbe1$hostID)
        } else {
          rewire_pullnodes_complete(pbe1$hostID)
        }
        propose_pbe("withinhost")
        if(runif(1) < pbe1$logLikseq - pbe0$logLikseq) {
          accept_pbe("withinhost")
        }
      }
    }
  }
  
  update_edge <- function(edgeid) {
    # tips in hostID
    host_of_edge <- pbe0$v$nodehosts[edgeid]
    
    # deconnect, reconnect, accept/reject 
    if(pbe0$v$nodeparents[edgeid] > 0) {
      prepare_pbe()
      coalnode <- take_cnode(edgeid)
      pbe1$v$nodehosts[edgeid] <- host_of_edge
      rewire_pullnodes_complete(pbe1$hostID)
      propose_pbe("withinhost")
      if(runif(1) < pbe1$logLikseq - pbe0$logLikseq) {
        accept_pbe("withinhost")
      }
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
    logaccprob <- pbe1$logLikgen + pbe1$logLiksam + pbe1$logLikcoal + pbe1$logLikdist - 
      pbe0$logLikgen - pbe0$logLiksam - pbe0$logLikcoal - pbe0$logLikdist + logproposalratio
    
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
    nodeUP <- .ptr(v$nodeparents, hostID-1)[v$nodetimes[.ptr(v$nodeparents, hostID-1)] < tinf.prop][1]
    nodeDP <- intersect(which(v$nodeparents == nodeUP), .ptr(v$nodeparents, hostID-1))
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
    v$infectors <- tail(v$nodehosts, p$obs+1)
    
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
