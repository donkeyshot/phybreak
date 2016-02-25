NumError <- 10^-14

###Key functions in updating process

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
      cnodes <- sort(nIDs[ntypes=="c"],decreasing = TRUE)
      res <- c(rep(NA,sum(ntypes=="s")),
               rootnode,tail(-nIDs[ntypes=="c"],-1),
               rep(NA,sum(ntypes=="t")))
      for(i in cnodes) {
        res[sample(which(is.na(res)),2)] <- i
        res[res == -i] <- NA
      }
      return(res)
    }
  )
  IDs <- nIDs[order(ntimes,ntypes)]
  tys <- ntypes[order(ntimes,ntypes)]
  if(tys[1] != "c") {
    print(c(nIDs, ntimes, ntypes, rootnode))
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



.lik.gentimes <- function(obs, shapeG, meanG, nodetimes, nodehosts, nodetypes) {
  sum(dgamma(nodetimes[nodetypes == "t" & nodehosts > 0] -
               nodetimes[nodehosts[nodetypes == "t" & nodehosts > 0] + 2*obs -1],
             shape = shapeG, scale = meanG/shapeG, log=TRUE))
}
.lik.sampletimes <- function(shapeS, meanS, nodetimes, nodetypes) {
  sum(dgamma(nodetimes[nodetypes == "s"] -
               nodetimes[nodetypes == "t"],
             shape = shapeS, scale = meanS/shapeS, log=TRUE))
}



###For initialization

#helper functions
.rinftimes <- function(st, meanS, shapeS) {
  st - rgamma(length(st), shape = shapeS, scale = meanS/shapeS)
}
.rinfectors <- function(it, meanG, shapeG) {
  if(sum(it == min(it)) > 1) stop("rinfectors with >1 index case")
  res <- rep(0,length(it))
  for(i in 1:length(it)) {
    if(it[i] > min(it)) {
      dist <- dgamma(it[i] - it, shape = shapeG, scale = meanG/shapeG)
      dist[i] <- 0
      res[i] <- sample(length(it), 1, prob = dist)
    }
  }
  return(res)
}
.distmatrix <- function(SNPs, SNPfreqs) {
  res <- matrix(0, nrow = nrow(SNPs), ncol = nrow(SNPs))
  for(i in 1:nrow(SNPs)) {
    for(j in i:nrow(SNPs)) {
      res[i,j] <- sum((SNPs[i,]!=SNPs[j,])*SNPfreqs)
      res[j,i] <- res[i,j]
    }
  }
  return(2^res)
}

#the actual initialization, given data, wh-model, and (initial) parameters, resulting in an phybreak-object
make.phybreak.obkData <- function(obk.object, mu = .01, shape.gen = 3, mean.gen = 1,
                                  shape.sample = 3, mean.sample = 1, wh.model = 3,
                                  wh.lambda= 0, wh.rate0 = 1,
                                  est.mean.gen = TRUE, est.mean.sample = TRUE,
                                  prior.mean.gen.shape = 0, prior.mean.gen.scale = 0,
                                  prior.mean.sample.shape = 0, prior.mean.sample.scale = 0,
                                  use.tree = FALSE) {
  obs <- nrow(obk.object@dna@dna[[1]])  #outbreaksize
  hostnames <- rownames(obk.object@dna@dna[[1]])

  refdate <- min(obk.object@dna@meta$date)

  if(!use.tree) {
    inftimes <- .rinftimes(as.numeric(obk.object@dna@meta$date - refdate), mean.sample, shape.sample)[1:obs]
    infectors <- .rinfectors(inftimes, mean.gen, shape.gen)
    nodeparents <- rep(0, 3*obs-1)
    nodetimes <- nodeparents   #initialize nodes: will contain time of node
    nodehosts <- nodeparents   #initialize nodes: will contain host carrying the node
    nodetypes <- c(rep("s",obs), rep("c",obs - 1),
                   rep("t",obs))  #initialize nodes: will contain node type (sampling, coalescent, transmission)
    nodetimes[1:obs] <- as.numeric(obk.object@dna@meta$date - refdate)[1:obs]   #sampling nodes at sample times
    nodetimes[1:obs + 2*obs - 1] <- inftimes  #transmission nodes at infection times
    nextcoalnode <- obs + 1  #needed in next loop: which is the next available coalescence node
    ## distribute nodes over the hosts
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
                         wh.model, wh.lambda, wh.rate0)  #...sampled coalescence times
    }
    ## sample for each node its parent node
    for(i in 1:obs) {
      nodeparents[nodehosts == i] <-     #change the parent nodes of all nodes in host i...
        .sampletopology(which(nodehosts == i), nodetimes[nodehosts == i], nodetypes[nodehosts == i], i + 2*obs - 1, wh.model)
      #...to a correct topology, randomized where possible
    }
    varlist <- list(
      nodetimes = nodetimes,
      nodehosts = nodehosts,
      nodeparents = nodeparents,
      nodetypes = nodetypes
    )
  } else {
    infectors <- obk.object@individuals$infector
    inftimes <- as.numeric(obk.object@individuals$date - refdate)
    samtimes <- as.numeric(obk.object@dna@meta$date - refdate)

    phylobegin <- obk.object@trees[[1]]$edge[,1]
    phyloend <- obk.object@trees[[1]]$edge[,2]
    phylolengths <- obk.object@trees[[1]]$edge.length

    nodeparents <- c(head(phylobegin[order(phyloend)],obs),
                     which(infectors==0)+2*obs-1,
                     tail(phylobegin[order(phyloend)],obs-2),
                     rep(NA,obs))
    nodetimes <- c(samtimes, rep(NA,obs-1), inftimes)   #initialize nodes: will contain time of node
    nodetimes <- c(rep(NA,2*obs-1), inftimes)   #initialize nodes: will contain time of node
    nodehosts <- c(1:obs, which(infectors==0), rep(NA,obs-2), infectors)   #initialize nodes: will contain host carrying the node
    edgelengths <- c(head(phylolengths[order(phyloend)],obs),
                     0, tail(phylolengths[order(phyloend)],obs-2),
                     rep(NA,obs))
    nodetypes <- c(rep("s",obs), rep("c",obs - 1),
                   rep("t",obs))  #initialize nodes: will contain node type (sampling, coalescent, transmission)
    while(any(is.na(nodetimes))) {
      nodetimes[1:(2*obs - 1)] <- nodetimes[nodeparents[1:(2*obs - 1)]] + edgelengths[1:(2*obs - 1)]
    }
    nodetimes <- round(nodetimes, digits = 12)
    edgelengths[1:(2*obs - 1)] <- nodetimes[1:(2*obs - 1)] - nodetimes[nodeparents[1:(2*obs - 1)]]
    
    nodeparents[(2*obs - 1) + setdiff(1:obs,nodehosts[(2*obs):(3*obs-1)])] <- nodeparents[setdiff(1:obs,nodehosts[(2*obs):(3*obs-1)])]
    nodeparents[setdiff(1:obs,nodehosts[(2*obs):(3*obs-1)])] <- (2*obs - 1) + setdiff(1:obs,nodehosts[(2*obs):(3*obs-1)])
    nodehosts[nodeparents[sort(unique(nodehosts[(2*obs):(3*obs-1)]))[-1]]] <- sort(unique(nodehosts[(2*obs):(3*obs-1)]))[-1]
    
    #nh[np[sameparentsamehost]] <- nh[sameparentsamehost]
    while(any(is.na(nodehosts))) {
      whichnodes <- tail(which(!is.na(nodehosts) & !is.na(nodeparents)),-obs)
      whichnodes <- whichnodes[order(nodeparents[whichnodes])]
      sameparent.wn <- which(duplicated(nodeparents[whichnodes]))
      for(i in sameparent.wn) {
        if(nodehosts[whichnodes[i]] == nodehosts[whichnodes[i-1]]) {
          nodehosts[nodeparents[whichnodes[i]]] <- nodehosts[whichnodes[i]]
        } else {
          posshosts <- nodehosts[whichnodes[i - 1:0]]
          posshosts <- c(posshosts, infectors[posshosts])
          nodehosts[nodeparents[whichnodes[i]]] <- posshosts[duplicated(posshosts)]
        }
      }
      
    }
    
    nodeparents[nodehosts[which(nodehosts[nodeparents[(obs+1):(2*obs-1)]] != nodehosts[(obs+1):(2*obs-1)])+obs]+2*obs-1] <- nodeparents[which(nodehosts[nodeparents[(obs+1):(2*obs-1)]] != nodehosts[(obs+1):(2*obs-1)])+obs]
    nodeparents[which(nodehosts[nodeparents[(obs+1):(2*obs-1)]] != nodehosts[(obs+1):(2*obs-1)])+obs] <- nodehosts[which(nodehosts[nodeparents[(obs+1):(2*obs-1)]] != nodehosts[(obs+1):(2*obs-1)])+obs]+2*obs-1
    nodeparents[nodehosts==0] <- 0
    

    varlist <- list(
      nodetimes = nodetimes,
      nodehosts = nodehosts,
      nodeparents = nodeparents,
      nodetypes = nodetypes
    )
  }

  SNP.sample <- c()
  SNP.frequencies <- c()
  for(i in 1:length(obk.object@dna@dna)) {
    SNP.sample <- cbind(SNP.sample,
                        do.call(rbind,as.phyDat(obk.object@dna@dna[[i]])))
    SNP.frequencies <- c(SNP.frequencies,
                         attr(as.phyDat(obk.object@dna@dna[[i]]),"weight"))
  }
  seq.sample <- SNP.sample
  for(i in (ncol(SNP.sample)-1):1) {
    for(j in length(SNP.frequencies):(i+1)) {
      if(all(seq.sample[,i] == seq.sample[,j])) {
        SNP.frequencies[i] <- SNP.frequencies[i] + SNP.frequencies[j]
        SNP.frequencies <- SNP.frequencies[-j]
        seq.sample <- seq.sample[,-j]
      }
    }
  }

  seq.sample[seq.sample == 1] <- "a"
  seq.sample[seq.sample == 2] <- "c"
  seq.sample[seq.sample == 3] <- "g"
  seq.sample[seq.sample == 4] <- "t"


  res <- list(
    d = list(
      names = hostnames,
      SNP = seq.sample,
      SNPfr = SNP.frequencies
    ),
    v = varlist,
    p = list(
      mu = mu,
      mean.gen = mean.gen,
      mean.sample = mean.sample,
      obs = obs,
      shape.gen = shape.gen,
      shape.sample = shape.sample,
      wh.model = wh.model,
      wh.lambda = wh.lambda,
      wh.rate0 = wh.rate0
    ),
    h = list(si = c(rep(c(mu,2*mu),50),rep(NA,400)),
             dist = .distmatrix(seq.sample, SNP.frequencies),
             est.mG = est.mean.gen,
             est.mS = est.mean.sample,
             mG.sh = prior.mean.gen.shape,
             mG.sc = prior.mean.gen.scale,
             mS.sh = prior.mean.sample.shape,
             mS.sc = prior.mean.sample.scale),
    s = list(
      nodetimes = c(),
      nodehosts = c(),
      nodeparents = c(),
      mu = c(),
      mean.gen = c(),
      mean.sample = c(),
      logLik = c()
    )
  )
  return(res)

}




make.phybreak <- function (t.sample, SNP.sample, SNP.freqs = c(), hostnames = NULL, mu = .01,
                           shape.gen = 3, mean.gen = 1, shape.sample = 3, mean.sample = 1,
                           wh.model = 3, wh.lambda= 0, wh.rate0 = 1,
                           est.mean.gen = TRUE, est.mean.sample = TRUE,
                           prior.mean.gen.shape = 0, prior.mean.gen.scale = 0,
                           prior.mean.sample.shape = 0, prior.mean.sample.scale = 0) {
  obs <- length(t.sample)  #outbreaksize
  if(length(hostnames) != obs) {
    if(is.null(names(t.sample))) {
      hostnames <- paste0("host.",1:obs)
    } else hostnames <- names(t.sample)
  }


  inftimes <- .rinftimes(t.sample,mean.sample,shape.sample)
  infectors <- .rinfectors(inftimes,mean.gen,shape.gen)
  nodeparents <- rep(0,3*obs - 1)  #initialize nodes: will contain parent node in phylotree
  nodetimes <- nodeparents   #initialize nodes: will contain time of node
  nodehosts <- nodeparents   #initialize nodes: will contain host carrying the node
  nodetypes <- c(rep("s",obs), rep("c",obs - 1),
                 rep("t",obs))  #initialize nodes: will contain node type (sampling, coalescent, transmission)
  nodetimes[1:obs] <- t.sample   #sampling nodes at sample times
  nodetimes[1:obs + 2*obs - 1] <- inftimes  #transmission nodes at infection times
  nextcoalnode <- obs + 1  #needed in next loop: which is the next available coalescence node
  ## distribute nodes over the hosts
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
                       wh.model, wh.lambda, wh.rate0)  #...sampled coalescence times
  }
  ## sample for each node its parent node
  for(i in 1:obs) {
    nodeparents[nodehosts == i] <-     #change the parent nodes of all nodes in host i...
      .sampletopology(which(nodehosts == i), nodetimes[nodehosts == i], nodetypes[nodehosts == i], i + 2*obs - 1, wh.model)
    #...to a correct topology, randomized where possible
  }
  varlist <- list(
    nodetimes = nodetimes,
    nodehosts = nodehosts,
    nodeparents = nodeparents,
    nodetypes = nodetypes
  )


  if(length(SNP.freqs) == 0) {
    SNP.frequencies <- rep(1,ncol(SNP.sample))
  } else {
    SNP.frequencies <- SNP.freqs
  }
  seq.sample <- SNP.sample
  for(i in (ncol(SNP.sample)-1):1) {
    for(j in length(SNP.frequencies):(i+1)) {
      if(all(seq.sample[,i] == seq.sample[,j])) {
        SNP.frequencies[i] <- SNP.frequencies[i] + SNP.frequencies[j]
        SNP.frequencies <- SNP.frequencies[-j]
        seq.sample <- seq.sample[,-j]
      }
    }
  }

  res <- list(
    d = list(
      names = hostnames,
      SNP = seq.sample,
      SNPfr = SNP.frequencies
    ),
    v = varlist,
    p = list(
      mu = mu,
      mean.gen = mean.gen,
      mean.sample = mean.sample,
      obs = obs,
      shape.gen = shape.gen,
      shape.sample = shape.sample,
      wh.model = wh.model,
      wh.lambda = wh.lambda,
      wh.rate0 = wh.rate0
    ),
    h = list(si = c(rep(c(mu,2*mu),50),rep(NA,400)),
             dist = .distmatrix(seq.sample, SNP.frequencies),
             est.mG = est.mean.gen,
             est.mS = est.mean.sample,
             mG.sh = prior.mean.gen.shape,
             mG.sc = prior.mean.gen.scale,
             mS.sh = prior.mean.sample.shape,
             mS.sc = prior.mean.sample.scale),
    s = list(
      nodetimes = c(),
      nodehosts = c(),
      nodeparents = c(),
      mu = c(),
      mean.gen = c(),
      mean.sample = c(),
      logLik = c()
    )
  )
  return(res)
}


clearsamples.phybreak <- function(phybreak.object, clearall = TRUE, nkeep = 1000) {
  if(clearall) {
    return(
      within(
        phybreak.object,
        s <- list(
          nodetimes = c(),
          nodehosts = c(),
          nodeparents = c(),
          mu = c(),
          mG = c(),
          mS = c(),
          logLik = c()
        )
      )
    )
  }
  if(length(phybreak.object$s$logLik) < nkeep) return(phybreak.object)
  tokeep <- length(phybreak.object$s$logLik) - (nkeep - 1):0
  return(
    within(
      phybreak.object,
      s <- list(
        nodetimes = s$nodetimes[,tokeep],
        nodehosts = s$nodehosts[,tokeep],
        nodeparents = s$nodeparents[,tokeep],
        mu = s$mu[tokeep],
        mG = s$mG[tokeep],
        mS = s$mS[tokeep],
        logLik = s$logLik[tokeep]
      )
    )
  )
}




