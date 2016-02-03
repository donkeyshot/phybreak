.samplecoaltimes <- function(tleaves, WHmodel = 3, lambda = 0, rate0 = 1) {
  if(length(tleaves) < 2) return(c())

  switch(
    WHmodel,
    return(head(sort(tleaves),-1)),
    return(0*tleaves[-1]),
    {
      # transform times so that fixed rate 1 can be used
      if(lambda == 0) {
        ttrans <- sort(tleaves/rate0,decreasing = TRUE)
      } else {
        ttrans <- sort((1-exp(-lambda*tleaves))/(lambda*rate0), decreasing = TRUE)
      }
      #       tnodetrans <- c(0)
      #       for( i in 1:length(ttrans)) {
      #         currentnodetime <- ttrans[i]    #starting at leaf i
      #         nodetimesinpast <- tnodetrans[tnodetrans < currentnodetime]
      #         totalrate <- sum(currentnodetime - nodetimesinpast)   #total coal rate to be exposed to
      #         cumratetocoalescence <- -log(runif(1,exp(-totalrate),1))          #conditional on coalescence within this host
      #         while(cumratetocoalescence > (currentnodetime - nodetimesinpast[1]) * length(nodetimesinpast)) {
      #           cumratetocoalescence <- cumratetocoalescence - (currentnodetime - nodetimesinpast[1]) * length(nodetimesinpast)
      #           currentnodetime <- nodetimesinpast[1]
      #           nodetimesinpast <- nodetimesinpast[-1]
      #         }
      #         tnodetrans <- sort(c(tnodetrans, currentnodetime - cumratetocoalescence / length(nodetimesinpast)), decreasing = TRUE)
      #       }
      tnodetrans <- .sct(ttrans)

      if(lambda == 0) {
        return(sort(rate0 * tnodetrans))
      } else {
        return(sort(-log(1-rate0*lambda*tnodetrans)/lambda))
      }
    }
  )
}

.sampletopology <- function(nIDs, ntimes, ntypes, rootnode, WHmodel = 3) {
  if(length(nIDs) == 1) return(rootnode)
  switch(
    WHmodel,{
      cnodes <- nIDs[ntypes=="c"]
      cnodeparents <- c(rootnode,head(cnodes,-1))
      leafparents <- c(cnodes,tail(cnodes,1))
      leafparents <- leafparents[rank(ntimes[ntypes != "c"],ties.method="first")]
      res <- c(head(leafparents,sum(ntypes=="s")),
               cnodeparents,
               tail(leafparents,sum(ntypes=="t")))
      return(res)
    },{
      cnodes <- nIDs[ntypes=="c"]
      cnodeparents <- c(rootnode,head(cnodes,-1))
      leafparents <- c(cnodes,tail(cnodes,1))
      res <- c(head(leafparents,sum(ntypes=="s")),
               cnodeparents,
               tail(leafparents,sum(ntypes=="t")))
      return(res)
    }
  )
  IDs <- nIDs[order(ntimes, ntypes)]
  tys <- ntypes[order(ntimes, ntypes)]
  if(tys[1] != "c") {
    stop("host topology does not start with coalescence node")
  }
  res <- rep(rootnode, length(nIDs))
  tochoose <- rep(IDs[1], 2)
  for(i in 2:length(nIDs)) {
    res[i] <- tochoose[1]
    if(tys[i] == "c") {
      tochoose <- sample(c(tochoose[-1], IDs[i], IDs[i]))
    } else {
      tochoose <- tochoose[-1]
    }
  }
  return(res[order(IDs)])
}


.sim.outbreak.gentime <- function(Npop = 100, R0 = 1.5, aG = 3, mG = 1, aS = 3, mS = 1,
                             nSUP = 0, multSUP = 10, disp = Inf) {
  inftimes <- c(0, rep(10000, Npop-1))

  infectivity <- if(disp== Inf) {
    sample(c(rep(multSUP * R0, nSUP), rep(R0, Npop - nSUP)))
  } else {
    sample(c(rgamma(nSUP,disp,disp/(multSUP*R0)), rgamma(Npop-nSUP,disp,disp/R0)))
  }
  nrcontacts <- rpois(Npop,infectivity)

  nth.infection <- 1
  currentID <- 1
  sources <- rep(0,Npop)

  while(nth.infection <= Npop & inftimes[currentID] != 10000) {
    whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
      #sorting so that with double contacts, the earliest will be used last
    whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE)
    successful <- whencontacts < inftimes[whocontacted]
    sources[whocontacted[successful]] <- currentID
    inftimes[whocontacted[successful]] <- whencontacts[successful]
    nth.infection <- nth.infection + 1
    currentID <- order(inftimes)[nth.infection]
  }

  obs <- sum(inftimes<10000)
  samptimes <- inftimes + rgamma(Npop, aS, aS/mS)

  orderbysamptimes <- order(samptimes)
  sources <- sources[orderbysamptimes]
  infectors <- match(sources,orderbysamptimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamptimes]
  samptimes <- samptimes[orderbysamptimes]

  return(
    list(
      obs = obs,
      sampletimes = samptimes[1:obs],
      infectiontimes = inftimes[1:obs],
      infectors = infectors
    )
  )
}

.sim.outbreak.gentime.size <- function(obsize = 50, R0 = 1.5, aG = 3, mG = 1, aS = 3, mS = 1,
                                      nSUP = 0, multSUP = 10, disp = Inf) {
  npop <- obsize
  while(1 - obsize/npop < exp(-(R0 * (1 + nSUP*(multSUP-1)/npop)) * obsize/npop)) {npop <- npop + 1}

  sim <- .sim.outbreak.gentime(npop, R0, aG, mG, aS, mS, nSUP, multSUP, disp)

  while(sim$obs != obsize) {
    sim <- .sim.outbreak.gentime(npop, R0, aG, mG, aS, mS, nSUP, multSUP, disp)
  }

  return(sim)
}

.sim.outbreak.SIR <- function(Npop = 100, R0 = 1.5, aI = 3, mI = 1, nSUP = 0, multSUP = 10) {
  inftimes <- c(0, rep(10000, Npop-1))

  infperiods <- sample(c(rgamma(nSUP,aI,aI/(mI*multSUP)), rgamma(Npop-nSUP,aI,aI/mI)))
  infectivity <- R0 * infperiods
  nrcontacts <- rpois(Npop,infectivity)

  nth.infection <- 1
  currentID <- 1
  sources <- rep(0,Npop)

  while(nth.infection <= Npop & inftimes[currentID] != 10000) {
    whencontacts <- sort(inftimes[currentID] + runif(nrcontacts[currentID], 0, infperiods[currentID]),decreasing = TRUE)
    #sorting so that with double contacts, the earliest will be used last
    whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE)
    successful <- whencontacts < inftimes[whocontacted]
    sources[whocontacted[successful]] <- currentID
    inftimes[whocontacted[successful]] <- whencontacts[successful]
    nth.infection <- nth.infection + 1
    currentID <- order(inftimes)[nth.infection]
  }

  obs <- sum(inftimes<10000)
  samptimes <- inftimes + infperiods

  orderbysamptimes <- order(samptimes)
  sources <- sources[orderbysamptimes]
  infectors <- match(sources,orderbysamptimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamptimes]
  samptimes <- samptimes[orderbysamptimes]

  return(
    list(
      obs = obs,
      sampletimes = samptimes[1:obs],
      infectiontimes = inftimes[1:obs],
      infectors = infectors
    )
  )

}

.sim.outbreak.SIR.size <- function(obsize = 50, R0 = 1.5, aI = 3, mI = 1, nSUP = 0, multSUP = 10) {
  npop <- obsize
  while(1 - obsize/npop < exp(-(R0 * (1 + nSUP*(multSUP-1)/obsize)) * obsize/npop)) {npop <- npop + 1}

  sim <- .sim.outbreak.SIR(npop, R0, aI, mI, nSUP, multSUP)

  while(sim$obs != obsize) {
    sim <- .sim.outbreak.SIR(npop, R0, aI, mI, nSUP, multSUP)
  }

  return(sim)
}


.sim.phylotree <- function (sim.object, wh.model = 3, lambda = 0, rate0 = 1) {
  with(
    sim.object,
    {
      nodeparents <- rep(0,3*obs - 1)  #initialize nodes: will containsparent node in phylotree
      nodetimes <- nodeparents   #initialize nodes: will contain time of node
      nodehosts <- nodeparents   #initialize nodes: will contain host carrying the node
      nodetypes <- c(rep("s",obs), rep("c",obs - 1),
                      rep("t",obs))  #initialize nodes: will contain node type (sampling, coalescent, transmission)

      nodetimes[1:obs] <- sampletimes   #sampling nodes at sample times
      nodetimes[1:obs + 2*obs - 1] <- infectiontimes  #transmission nodes at infection times
      nextcoalnode <- obs + 1  #needed in next loop: which is the next available coalescence node

      ## distribute nodes over the hosts
      nodehosts[1:obs] <- 1:obs
      nodehosts[(obs+1):(2*obs-1)] <- sort(infectors)[-1]
      nodehosts[(2*obs):(3*obs-1)] <- infectors
      for(i in 1:obs) {
        iees <- which(i == infectors)  #infectees of host i (transmission nodes in host i)
        nodehosts[i] <- i    #assign sample node
        while(length(iees) > 0) {   #per infectee...
          nodehosts[iees[1] + 2*obs - 1] <- i   #...assign transmission node
          nodehosts[nextcoalnode] <- i         #...assign coalescence node
          nextcoalnode <- nextcoalnode + 1     #next coalescence node
          iees <- iees[-1]                    #infectee done with
        }
      }

      ## sample the times of the coalescence nodes
      for(i in 1:obs) {
        nodetimes[nodehosts == i & nodetypes == "c"] <-   #change the times of the coalescence nodes in host i...
          nodetimes[i + 2*obs - 1] +                      #...to the infection time +
          .samplecoaltimes(nodetimes[nodehosts == i & nodetypes != "c"] - nodetimes[i + 2*obs - 1],
                          wh.model, lambda, rate0)  #...sampled coalescence times
      }
      ## sample for each node its parent node
      for(i in 1:obs) {
        nodeparents[nodehosts == i] <-     #change the parent nodes of all nodes in host i...
          .sampletopology(which(nodehosts == i), nodetimes[nodehosts == i], nodetypes[nodehosts == i], i + 2*obs - 1, wh.model)
        #...to a correct topology, randomized where possible
      }
      return(
        within(sim.object, {
          nodetypes <- nodetypes
          nodeparents <- nodeparents
          nodehosts <- nodehosts
          nodetimes <- nodetimes
        }))


    })

}

.sim.sequences <- function (sim.object, mutrate, nr.sites = 10000) {
  with(sim.object,{
    # sim.object,{
    edgelengths <- nodetimes - c(0,nodetimes)[1+nodeparents]
    edgelengths[edgelengths < 0] <- 0
    nmutations <- rpois(1, mutrate*sum(edgelengths))
    mutedges <- sample(3*obs-1, size = nmutations, replace = TRUE, prob = edgelengths)
    mutedges <- mutedges[order(nodetimes[mutedges])]
    mutsites <- sample(nr.sites, size = nmutations, replace = TRUE)
    mutsites <- match(mutsites, unique(mutsites))
    mutnucl <- sample(4, size = nmutations, replace = TRUE)


    nodestrains <- matrix(data = rep(1, nr.sites * (3*obs-1)), ncol = nr.sites)
    for(i in 1:(3*obs-1)) {
      currentedge <- i
      while(nodeparents[currentedge] != 0) {
        nodestrains[i,mutsites[mutedges == currentedge]] <-
          mutnucl[mutedges == currentedge]
        currentedge <- nodeparents[currentedge]
      }

    }
    nodestrains <- nodestrains[nodetypes == "s",]
    #       if(length(unique(mutsites)) < nr.sites) {
    #         nodestrains <- nodestrains[, 1:(length(unique(mutsites))+1)]
    #       }

    #       SNPfrequencies <- c(rep(1,length(unique(mutsites))), nr.sites - length(unique(mutsites)))
    #       if(tail(SNPfrequencies,1) == 0) SNPfrequencies <- head(SNPfrequencies, -1)
    #       for(i in (length(SNPfrequencies)-1):1) {
    #         for(j in length(SNPfrequencies):(i+1)) {
    #           if(all(nodestrains[,i] == nodestrains[,j])) {
    #             SNPfrequencies[i] <- SNPfrequencies[i] + SNPfrequencies[j]
    #             SNPfrequencies <- SNPfrequencies[-j]
    #             nodestrains <- nodestrains[,-j]
    #           }
    #         }
    #       }

    nodestrains[nodestrains == 1] <- "a"
    nodestrains[nodestrains == 2] <- "c"
    nodestrains[nodestrains == 3] <- "g"
    nodestrains[nodestrains == 4] <- "t"

    rownames(nodestrains) <- 1:obs

    return(
      within(sim.object,{
        SNPlist <- as.DNAbin(nodestrains)
        #           SNP.frequencies <- SNPfrequencies
      })
    )
  }
  )
}

.makephylo <- function(nodetimes, nodeparents) {
  ###topology
  Nhosts <- (1+length(nodetimes))/3
  indext <- (1:length(nodetimes))[nodeparents == 0]
  indexc <- (1:length(nodetimes))[nodeparents == indext]
  edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
  edgeend <- (1:length(nodetimes))[nodeparents != 0 & nodeparents != indext]
  edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]

  toremove <- (edgeend >= 2*Nhosts)
  toreplace <- (edgestart >= 2*Nhosts)
  edgelengths[toreplace] <- edgelengths[toreplace] + edgelengths[toremove][match(edgestart[toreplace],edgeend[toremove])]
  edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace],edgeend[toremove])]
  edgestart <- edgestart[!toremove]
  edgeend <- edgeend[!toremove]
  edgelengths <- edgelengths[!toremove]
  if(indexc != Nhosts + 1) {
    edgestart[edgestart == indexc] <- 0
    edgeend[edgeend == indexc] <- 0
    edgestart[edgestart == Nhosts + 1] <- indexc
    edgeend[edgeend == Nhosts + 1] <- indexc
    edgestart[edgestart == 0] <- Nhosts + 1
    edgeend[edgeend == 0] <- Nhosts + 1
  }

  edges <- matrix(c(edgestart,edgeend),ncol=2)



  res <- list(
    edge = edges,
    edge.length = edgelengths,
    Nnode = Nhosts - 1,
    tip.label = 1:Nhosts
  )
  class(res) <- "phylo"
  res <- reorder(res)
  return(res)

}


sim.phybreak.gentime <- function(obsize = 50, R0 = 1.5, shape.gen = 3,
                                 mean.gen = 1, shape.sample = 3, mean.sample = 1,
                                 nSUP = 0, multSUP = 10, disp = Inf, wh.model = 3,
                                 lambda = 0, rate0 = 1, mutrate = 1, nr.sites = 10000) {
  res <- .sim.outbreak.gentime.size(obsize,R0,shape.gen,mean.gen,
                                    shape.sample,mean.sample,nSUP,multSUP,disp)
  res <- .sim.phylotree(res,wh.model,lambda,rate0)
  res <- .sim.sequences(res,mutrate,nr.sites)

  treesout <- vector('list',1)
  treesout[[1]] <- .makephylo(res$nodetimes, res$nodeparents)
  class(treesout) <- "multiPhylo"

  toreturn <- new("obkData",
                  individuals = data.frame(
                    infector = res$infectors,
                    date = as.Date(res$infectiontimes, origin = "2000-01-01"),
                    row.names = 1:res$obs),
                  dna = list(SNPs = res$SNPlist), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                  dna.individualID = 1:res$obs, trees = treesout)
  return(toreturn)
#   return(list(
#     dataset = list(
#       obs = res$obs,
#       sampletimes = res$sampletimes,
#       SNPs = res$SNPlist#,
# #       SNP.frequencies = res$SNP.frequencies
#     ),
#     trueoutbreak = list(
#       nodetimes = res$nodetimes,
#       nodehosts = res$nodehosts,
#       nodeparents = res$nodeparents,
#       nodetypes = res$nodetypes
#     )
#   ))
}

sim.phybreak.SIR <- function(obsize = 50, R0 = 1.5, shape.infper = 3, mean.infper = 1,
                                 nSUP = 0, multSUP = 10, wh.model = 3, lambda = 0, rate0 = 1,
                                 mutrate = 1, nr.sites = 10000) {
  res <- .sim.outbreak.SIR.size(obsize,R0,shape.infper,mean.infper,nSUP,multSUP)
  res <- .sim.phylotree(res,wh.model,lambda,rate0)
  res <- .sim.sequences(res,mutrate,nr.sites)

  treesout <- vector('list',1)
  treesout[[1]] <- .makephylo(res$nodetimes, res$nodeparents)
  class(treesout) <- "multiPhylo"

  toreturn <- new("obkData",
                  individuals = data.frame(
                    infector = res$infectors,
                    date = as.Date(res$infectiontimes, origin = "2000-01-01"),
                    row.names = 1:res$obs),
                  dna = list(SNPs = res$SNPlist), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                  dna.individualID = 1:res$obs, trees = treesout)
  return(toreturn)
  #   return(list(
  #     dataset = list(
  #       obs = res$obs,
  #       sampletimes = res$sampletimes,
  #       SNPs = res$SNPlist#,
  # #       SNP.frequencies = res$SNP.frequencies
  #     ),
  #     trueoutbreak = list(
  #       nodetimes = res$nodetimes,
  #       nodehosts = res$nodehosts,
  #       nodeparents = res$nodeparents,
  #       nodetypes = res$nodetypes
  #     )
  #   ))
}

resim.phybreak <- function(obk.object, wh.resim = TRUE, seq.resim = TRUE,
                           wh.model = 3, lambda = 0, rate0 = 1,
                           mutrate = 1, nr.sites = 10000) {
  refdate <- min(obk.object@individuals$date)
  if(wh.resim) {
    res <- list(
      obs = length(obk.object@individuals$infector),
      sampletimes = as.numeric(obk.object@dna@meta$date - refdate),
      infectiontimes = as.numeric(obk.object@individuals$date - refdate),
      infectors = obk.object@individuals$infector
    )
    res <- .sim.phylotree(res,wh.model,lambda,rate0)
    res <- .sim.sequences(res,mutrate,nr.sites)

    treesout <- vector('list',1)
    treesout[[1]] <- .makephylo(res$nodetimes, res$nodeparents)
    class(treesout) <- "multiPhylo"

    toreturn <- new("obkData",
                    individuals = data.frame(
                      infector = res$infectors,
                      date = as.Date(res$infectiontimes, origin = "2000-01-01"),
                      row.names = 1:res$obs),
                    dna = list(SNPs = res$SNPlist), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                    dna.individualID = 1:res$obs, trees = treesout)
    return(toreturn)
  } else if (seq.resim) {
    tempres <- make.phybreak.obkData(obk.object)

    res <- list(
      obs = length(obk.object@individuals$infector),
      sampletimes = as.numeric(obk.object@dna@meta$date - refdate),
      infectiontimes = as.numeric(obk.object@individuals$date - refdate),
      infectors = obk.object@individuals$infector,
      nodetypes = tempres$v$nodetypes,
      nodeparents = tempres$v$nodeparents,
      nodehosts = tempres$v$nodehosts,
      nodetimes = tempres$v$nodetimes - min(tempres$v$nodetimes)
    )
    res <- .sim.sequences(res,mutrate,nr.sites)

    treesout <- vector('list',1)
    treesout[[1]] <- .makephylo(res$nodetimes, res$nodeparents)
    class(treesout) <- "multiPhylo"

    toreturn <- new("obkData",
                    individuals = data.frame(
                      infector = res$infectors,
                      date = as.Date(res$infectiontimes, origin = "2000-01-01"),
                      row.names = 1:res$obs),
                    dna = list(SNPs = res$SNPlist), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                    dna.individualID = 1:res$obs, trees = treesout)
    return(toreturn)

  } else return(obk.object)

}


