#' Outbreak simulation.
#' 
#' Simulate outbreaks of class \code{phybreakdata}, with the outbreak model of \pkg{phybreak} (\code{sim.phybreak} is deprecated).
#' 
#' @param obsize The outbreak size (number of cases) to obtain. If \code{obsize = NA}, \code{popsize} should be provided.
#' @param popsize The population size in which to simulate. If it is not defined (default), 
#'   an optimal population size will be chosen based on R0 and obsize. Be aware that choosing a \code{popsize} and
#'   an \code{obsize} can severely increase the simulation time, depending on \code{R0}.
#' @param samplesperhost Number of samples to be taken per host, either a vector or a single number.
#' @param R0 The basic reproduction ratio used for simulation. The offspring distribution is Poisson.
#' @param introductions The number of index cases when simulating a multiple introduced outbreak.
#' @param introdistribution The probability distribution for the infection times of the index cases. Options are 
#'   \code{introdistribution = "exp"} for an exponential grow over time of the probability density, or 
#'   \code{introdistribution = "unif"} for an uniform probability distribution.
#' @param spatial If \code{TRUE}, the hosts are placed on a square with density 1, and a distance kernel is used to 
#'   model transmission probabilities between the hosts.
#' @param gen.shape The shape parameter of the gamma-distributed generation interval.
#' @param gen.mean The mean generation interval.
#' @param sample.shape The shape parameter of the gamma-distributed sampling interval.
#' @param sample.mean The mean sampling interval (for the first sample of each host).
#' @param additionalsampledelay Sampling intervals since first sampling times of each host. Values in this vector will be 
#'   used first for all additional samples of host 1, then of host 2, etc.
#' @param wh.model The model for within-host pathogen dynamics (effective pathogen population size = 
#'   N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Names and numbers are allowed.
#'   Options are:
#'   \enumerate{
#'     \item "single": effective size = 0, so coalescence occurs 'just before' transmission in the infector (complete bottleneck)
#'     \item "infinite": effective size = Inf, with complete bottleneck, so coalescence occurs 'just after' transmission in the infectee
#'     \item "linear": effective size at time t after infection = \code{wh.level + wh.slope * t} (complete or wide bottleneck; if complete, \code{wh.level = 0})
#'     \item "exponential": effective size at time t after infection = \code{wh.level * exp(wh.exponent * t)} (wide bottleneck)
#'     \item "constant": effective size = wh.level (wide bottleneck)
#'   }
#' @param wh.bottleneck Whether the bottleneck should be complete or wide, which is only an option if \code{wh.model = "linear"} 
#'   (in that case, \code{"auto"} defaults to \code{"complete"}).
#' @param wh.slope Within-host slope, used if \code{wh.model = "linear"}.
#' @param wh.exponent Within-host exponent, used if \code{wh.model = "exponential"}
#' @param wh.level Within-host effective pathogen size at transmission, used if \code{wh.bottleneck = "wide"}
#'   (if \code{wh.model = "exponential"} or \code{"constant"}, and optional if \code{wh.model = "linear"})
#' @param dist.model The distance kernel to use if \code{spatial = TRUE}. Options are:
#'   \enumerate{
#'     \item "power": a power law function pr(dist) ~ 1 / (1 + (dist/dist.scale) ^ dist.exponent)
#'     \item "exponential": an exponential function pr(dist) ~ exp(-dist.exponent * dist)
#'   }
#' @param dist.exponent Distance model exponent.
#' @param dist.scale Distance model scale, only with power law distance model.
#' @param mu Expected number of mutations per nucleotide per unit of time along each lineage. 
#' @param sequence.length Number of available nucleotides for mutations.
#' @param ... If arguments from previous versions of this function are used, they may be interpreted correctly through 
#'   this argument, but it is better to provide the correct argument names.
#' @return The simulation output as an object of class \code{'phybreakdata'} with sequences (class \code{'phyDat'}) and 
#'   sampling times (which would be the observations), and infection times, infectors, and phylogenetic tree 
#'   of class \code{\link[ape]{phylo}}.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' simulation <- sim_phybreak()
#' @export
sim_phybreak <- function(obsize = 50, popsize = NA, samplesperhost = 1,
                         introductions = 1, introdistribution = "unif", outbreaktime = 10,
                         R0 = 1.5, spatial = FALSE,
                         gen.shape = 10, gen.mean = 1, 
                         sample.shape = 10, sample.mean = 1, histtime = -10,
                         additionalsampledelay = 0,
                         wh.model = "linear", wh.bottleneck = "auto", wh.slope = 1, wh.exponent = 1, wh.level = 0.1,
                         dist.model = "power", dist.exponent = 2, dist.scale = 1,
                         mu = 0.0001, sequence.length = 10000, ...) {
  ### parameter name compatibility 
  old_arguments <- list(...)
  if(exists("old_arguments$shape.gen")) gen.shape <- old_arguments$shape.gen
  if(exists("old_arguments$mean.gen")) gen.mean <- old_arguments$mean.gen
  if(exists("old_arguments$shape.sample")) sample.shape <- old_arguments$shape.sample
  if(exists("old_arguments$mean.sample")) sample.mean <- old_arguments$sample.mean

  
  ### tests
  if(all(is.na(c(obsize, popsize)))) {
    stop("give an outbreak size (obsize) and/or a population size (popsize)")
  }
  if(all(!is.na(c(obsize, popsize)))) {
    warning("giving both an outbreak size (obsize) and a population size (popsize) can take a very large simulation time",
            immediate. = TRUE)
  }
  if(all(!is.na(c(obsize, popsize))) && obsize > popsize) {
    stop("outbreak size (obsize) cannot be larger than population size (popsize)")
  }
  if(R0 <= 1) {
    stop("R0 should be larger than 1")
  }
  if(any(c(gen.shape, gen.mean, sample.shape, sample.mean, wh.slope, mu) <= 0)) {
    stop("parameter values should be positive")
  }
  if(introdistribution != "exp" & introdistribution != "unif"){
    stop('introdistribution should be either "exp" or "unif"')
  }
  if(introductions == 1){
    histtime <- NULL
  }
  wh.model <- choose_whmodel(wh.model)
  wh.bottleneck <- choose_whbottleneck(wh.bottleneck, wh.model)
  wh.level <- wh.level * (wh.bottleneck == "wide")
  
  ### simulate step by step
  if(is.na(obsize)) {
    if(spatial) {
      res <- sim_outbreak_spatial(popsize, R0, introductions, introdistribution, gen.shape, gen.mean,
                          sample.shape, sample.mean, dist.model, dist.exponent, dist.scale)
    } else {
      res <- sim_outbreak(popsize, R0, introductions, introdistribution, gen.shape, gen.mean,
                          sample.shape, sample.mean)
    }
     obsize <- res$obs
     if(obsize == 1) return(c("Outbreak size = 1"))
  } else {
    if(spatial) {
      res <- sim_outbreak_size_spatial(obsize, popsize, R0, introductions, introdistribution, gen.shape, gen.mean,
                               sample.shape, sample.mean, dist.model, dist.exponent, dist.scale)
    } else {
      res <- sim_outbreak_size(obsize, popsize, R0, introductions, introdistribution, outbreaktime, gen.shape, gen.mean,
                               sample.shape, sample.mean)
    }
  }
  if(any(samplesperhost < 1)) stop("samplesperhost should be positive")
  if(any(additionalsampledelay < 0)) stop("additionalsampledelay cannot be negative")
  res <- sim_additionalsamples(res, samplesperhost, additionalsampledelay)
  res <- sim_phylotree(res, wh.model, wh.bottleneck, wh.slope, wh.exponent, wh.level, sample.mean, histtime)
  res <- sim_sequences(res, mu, sequence.length)
  
  hostnames <- paste0("host.", 1:obsize)
  
  if (introductions > 1) {
    samplenames <- paste0("sample.", res$nodehosts[1:res$Nsamples]-1, ".", nthsample(res))
  } else { 
    samplenames <- paste0("sample.", res$nodehosts[1:res$Nsamples], ".", nthsample(res))
  }
  names(res$sequences) <- samplenames
  if(spatial) {
    rownames(res$locations) <- hostnames
  }

  ### make a phylo tree
  treeout <- phybreak2phylo(vars = environment2phybreak(res), samplenames = samplenames, simmap = FALSE)
  
  if(introductions > 1){
    res <- within(res, {
      inftimes <- inftimes[-1]
      infectors <- infectors[-1]-1
    })
  }
  
  if(spatial) {
    toreturn <- with(res,
                     phybreakdata(
                       sequences = sequences,
                       sample.times = c(samtimes, addsampletimes),
                       spatial = locations,
                       sample.names = samplenames,
                       host.names = hostnames[c(1:obs, addsamplehosts)],
                       sim.infection.times = inftimes,
                       sim.infectors = infectors,
                       sim.tree = treeout,
                       external.seq = F
                ))
  } else {
    toreturn <- with(res,
                     phybreakdata(
      sequences = sequences,
      sample.times = c(samtimes, addsampletimes),
      sample.names = samplenames,
      host.names = hostnames[c(1:obs, addsamplehosts)],
      sim.infection.times = inftimes,
      sim.infectors = infectors,
      sim.tree = treeout,
      external.seq = F
    ))
  }
  
  return(toreturn)
}

#' @rdname sim_phybreak
#' @export
sim.phybreak <- function(...) {
  .Deprecated("sim_phybreak")
  sim_phybreak(...)
}


### simulate an outbreak of a particular size by repeating
### simulations until obsize is obtained
sim_outbreak_size <- function(obsize, Npop, R0, intro, introdist, time, aG, mG, aS, mS) {
  if(is.na(Npop)) {
    Npop <- obsize
    while(1 - obsize/Npop < exp(-R0* obsize/Npop)) {Npop <- Npop + 1}
  } 
  
  sim <- sim_outbreak(Npop, R0, intro, introdist, time, aG, mG, aS, mS)
  
  while(sim$obs != obsize) {
    sim <- sim_outbreak(Npop, R0, intro, introdist, time, aG, mG, aS, mS)
  }
  
  return(sim)
}

sim_outbreak_size_spatial <- function(obsize, Npop, R0, intro, introdist, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale) {
  if(is.na(Npop)) {
    Npop <- obsize
    while(1 - obsize/Npop < exp(-R0* obsize/Npop)) {Npop <- Npop + 1}
  } 
  
  sim <- sim_outbreak_spatial(Npop, R0, intro, introdist, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale)
  
  while(sim$obs != obsize | sum(sim$infectors==1) != intro) {
    sim <- sim_outbreak_spatial(Npop, R0, intro, introdist, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale)
  }
  
  return(sim)
}


### simulate an outbreak
sim_outbreak <- function(Npop, R0, intro, introdist, time, aG, mG, aS, mS) {
  ### initialize
  introwaittimes <- rexp(intro-1, rate = intro/time)
  introtimes <- c(0,cumsum(introwaittimes))
  
  Npops <- rep(1, intro)
  Npops_plus <- sample(intro, Npop - intro, replace = TRUE)
  Npops[sort(unique(Npops_plus))] <- Npops[sort(unique(Npops_plus))] + table(Npops_plus)
  
  trees <- lapply(1:intro, function(i){
    inftimes <- c(0, rep(10000, Npops[i]-1))
    sources <- rep(0,Npops[i])
    nrcontacts <- rpois(Npops[i], R0)
    nth.infection <- 1
    currentID <- 1
    
    # if (introdist == "exp")
    #   r <- (R0^(1/aG)-1)/(mG/aG)
    # else 
    #   r <- 0
    
    ### by order of infection, sample secondary infections
    # currentID is the infected host under consideration
    # while(nth.infection <= Npops[i] & inftimes[currentID] != 10000) {
    #   #when does currentID make infectious contacts?
    #   #reverse sorting so that with double contacts, the earliest will be used last
    #   # if (intro > 1){
    #   #   if (currentID == 1){
    #   #     #whencontacts <- sort(runif(nrcontacts[currentID],inftimes[currentID],inftimes[currentID]+log(Npop,R0)*mG),decreasing = TRUE)
    #   #     X <- runif(nrcontacts[currentID]-1)
    #   #     Y <- ifelse(rep(r,length(X))>0, log(X*exp(r*log(Npop,R0)*mG)-X+1)/r, runif(nrcontacts[currentID]-1,0,10))
    #   #     whencontacts <- sort(c(0,Y), decreasing = TRUE)
    #   #   } else {
    #   #     whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
    #   #   }
    #   # } else {
    #     whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
    #   # }
    #   
    #   #who are these contacts made with?
    #   whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE)
    #   
    #   #are these contacts successful, i.e. earlier than the existing contacts with these hosts?
    #   successful <- whencontacts < inftimes[whocontacted]
    #   
    #   #change infectors and infection times of successful contactees
    #   sources[whocontacted[successful]] <- currentID
    #   inftimes[whocontacted[successful]] <- whencontacts[successful]
    #   
    #   #go to next infected host in line
    #   nth.infection <- nth.infection + 1
    #   currentID <- order(inftimes)[nth.infection]
    # }
    # 
    while(nth.infection <= Npops[i] & inftimes[currentID] != 10000) {
      #when does currentID make infectious contacts?
      #reverse sorting so that with double contacts, the earliest will be used last
      whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
      
      #who are these contacts made with?
      whocontacted <- sample(Npops[i], nrcontacts[currentID], replace = TRUE)
      
      #are these contacts successful, i.e. earlier than the existing contacts with these hosts?
      successful <- whencontacts < inftimes[whocontacted]
      
      #change infectors and infection times of successful contactees
      sources[whocontacted[successful]] <- currentID
      inftimes[whocontacted[successful]] <- whencontacts[successful]
      
      #go to next infected host in line
      nth.infection <- nth.infection + 1
      currentID <- order(inftimes)[nth.infection]
    }
    
    ### determine outbreaksize and sampling times
    obs <- sum(inftimes<10000)
    # if (intro > 1)
    #   samtimes <- c(0,inftimes[-1] + rgamma(Npop, aS, aS/mS))
    # else 
    samtimes <- inftimes + rgamma(Npops[i], aS, aS/mS)
    
    ### order hosts by sampling times, and renumber hostIDs
    ### so that the uninfected hosts can be discarded
    orderbysamtimes <- order(samtimes)
    sources <- sources[orderbysamtimes]
    infectors <- match(sources,orderbysamtimes)[1:obs]
    infectors[is.na(infectors)] <- 0
    inftimes <- inftimes[orderbysamtimes]
    samtimes <- samtimes[orderbysamtimes]
    
    ### return the outbreak
    # if(intro > 1){
    #   inftimes <- inftimes - min(inftimes[-1])
    #   samtimes <- samtimes + inftimes[1]
    #   
    #   return(
    #     list(
    #       obs = obs-1,
    #       samtimes = samtimes[2:obs],
    #       inftimes = inftimes[2:obs],
    #       infectors = infectors
    #     )
    #   )
    # } else {
      return(
        list(
          ntree = i,
          obs = obs,
          samtimes = samtimes[1:obs] + introtimes[i],
          inftimes = inftimes[1:obs] + introtimes[i],
          infectors = infectors
        )
      )
  })
  
  hostorder <- order(do.call(c,lapply(trees, function(t) return(t$samtimes))))
  samtimes = do.call(c,lapply(trees, function(t) return(t$samtimes)))[hostorder]
  inftimes = do.call(c,lapply(trees, function(t) return(t$inftimes)))[hostorder]
  
  hosts <- match(1:length(hostorder), hostorder)
  
  infectors <- c()
  for (i in 1:intro){
    if (trees[[i]]$obs > 1){
      treehost <- hosts[1:trees[[i]]$obs]
      infectors <- c(infectors, 0, treehost[trees[[i]]$infectors[-1]])
    } else {
      infectors <- c(infectors, 0)
    }
    hosts <- hosts[-(1:trees[[i]]$obs)]
  }
  infectors <- infectors[hostorder]
  
  return(
    list(
      obs = length(infectors),
      samtimes = samtimes,
      inftimes = inftimes,
      infectors = c(0,infectors+1)
    )
  )
}

### simulate a spatial outbreak
sim_outbreak_spatial <- function(Npop, R0, intro, aG, mG, aS, mS, dist.model, dist.exponent, dist.scale) {
  ### initialize spatial population model
  x <- runif(Npop, 0, sqrt(Npop))
  y <- runif(Npop, 0, sqrt(Npop))
  distances <- as.matrix(dist(cbind(x, y)))
  dist_densities <- if(dist.model == "exponential") {
    dist.exponent * exp(-dist.exponent * distances)
  } else {
    dist.exponent * sin(pi/dist.exponent) / (pi * dist.scale * (1 + (distances/dist.scale) ^ dist.exponent))
  }
  dist_densities[cbind(1:Npop, 1:Npop)] <- 0
  # matrixR0 <- max(eigen(dist_densities)$values)
  matrixR0 <- sum(dist_densities)/Npop
  R0_matrix <- R0 * dist_densities / matrixR0
  if(intro > 1){
    R0_matrix <- cbind(0,R0_matrix)
    R0_matrix <- rbind(0,R0_matrix)
  }

  ### initialize outbreak
  inftimes <- c(0, rep(10000, Npop))
  sources <- rep(0, Npop+1)
  nrcontacts <- c(intro,rpois(Npop, rowSums(R0_matrix)))
  nth.infection <- 1
  currentID <- 1
  
  if (introdist == "exp")
    r <- (R0^(1/aG)-1)/(mG/aG)
  else 
    r <- 0
  
  ### by order of infection, sample secondary infections
  # currentID is the infected host under consideration
  while(nth.infection <= Npop+1 & inftimes[currentID] != 10000) {
    #when does currentID make infectious contacts?
    #reverse sorting so that with double contacts, the earliest will be used last
    whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
    
    #who are these contacts made with?
    if(intro > 1){
      if(currentID == 1) {
        whocontacted <- sample(Npop+1, nrcontacts[currentID], replace = FALSE)
        X <- runif(nrcontacts[currentID]-1)
        Y <- ifelse(rep(R0,length(X))>1, log(X*exp(r*log(Npop,R0)*mG)-X+1)/r, runif(nrcontacts[currentID]-1,0,log(Npop,R0)*mG))
        whencontacts <- sort(c(0,Y), decreasing = TRUE)
      } else {
        whocontacted <- sample(Npop+1, nrcontacts[currentID], replace = TRUE, prob = R0_matrix[currentID, ])
      }
    } else {
      whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE, prob = R0_matrix[currentID, ])
    }
    
    #are these contacts successful, i.e. earlier than the existing contacts with these hosts?
    successful <- whencontacts < inftimes[whocontacted]
    
    #change infectors and infection times of successful contactees
    sources[whocontacted[successful]] <- currentID
    inftimes[whocontacted[successful]] <- whencontacts[successful]
    #go to next infected host in line
    nth.infection <- nth.infection + 1
    currentID <- order(inftimes)[nth.infection]
  }
  
  ### determine outbreaksize and sampling times
  obs <- sum(inftimes<10000)
  if (intro > 1)
    samtimes <- c(0,inftimes[-1] + rgamma(Npop, aS, aS/mS))
  else 
    samtimes <- inftimes + rgamma(Npop+1, aS, aS/mS)
  
  ### order hosts by sampling times, and renumber hostIDs
  ### so that the uninfected hosts can be discarded
  orderbysamtimes <- order(samtimes)
  sources <- sources[orderbysamtimes]
  infectors <- match(sources,orderbysamtimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamtimes]
  samtimes <- samtimes[orderbysamtimes]
  locations <- cbind(x, y)[orderbysamtimes[-1]-1, ]

  ### return the outbreak
  if(intro > 1){
    inftimes <- inftimes - min(inftimes[-1])
    samtimes <- samtimes + inftimes[1]
    
    return(
      list(
        obs = obs-1,
        samtimes = samtimes[2:obs],
        inftimes = inftimes[2:obs],
        infectors = infectors,
        locations = locations[2:obs, ]
      )
    )
  } else {
    return(
      list(
        obs = obs,
        samtimes = samtimes[1:obs],
        inftimes = inftimes[1:obs],
        infectors = infectors,
        locations = locations[1:obs, ]
      )
    )
  }
}


### simulate additional samples in a transmission tree
sim_additionalsamples <- function(sim.object, samperh, addsamdelay) {
  # recycle too short arguments
  addsamplesizes <- rep_len(samperh - 1, sim.object$obs)
  addsamplesizes[addsamplesizes < 0] <- 0
  alldelays <- rep_len(addsamdelay, sum(addsamplesizes))
  
  # vectors with additional samplehosts and sample times
  addsamplehosts <- rep(1:sim.object$obs, addsamplesizes)
  addsampletimes <- sim.object$samtimes[addsamplehosts] + alldelays
  addsampletimes <- addsampletimes[order(addsamplehosts, addsampletimes)]
  
  return(within(sim.object, {
    Nsamples <- sim.object$obs + sum(addsamplesizes)
    addsamplehosts <- addsamplehosts
    addsampletimes <- addsampletimes
  }))
}

### simulate a phylogenetic tree given a transmission tree
sim_phylotree <- function (sim.object, wh.model, wh.bottleneck, wh.slope, wh.exponent, wh.level, sample.mean, histtime) {
  sim.object$inftimes <- c(histtime, sim.object$inftimes)
  list2env(list(v = sim.object, 
                p = list(wh.model = wh.model, wh.bottleneck = wh.bottleneck, wh.slope = wh.slope, wh.exponent = wh.exponent,
                         wh.level = wh.level, sample.mean = sample.mean),
                d = list(nsamples = sim.object$Nsamples)), pbe1)
  
  if(length(sim.object$infectors) != sim.object$obs){
    pbe1$v$nodeparents <- rep(-1, 2 * sim.object$Nsamples + sim.object$obs)  #initialize nodes: will contain parent node in phylotree
    pbe1$v$nodetimes <- c(sim.object$samtimes, sim.object$addsampletimes, 
                          rep(0, sim.object$Nsamples - 1), sim.object$inftimes)   #initialize nodes: will contain time of node
    pbe1$v$nodehosts <- c(2:(sim.object$obs+1), sim.object$addsamplehosts, 
                          rep(-1, sim.object$Nsamples - 1), sim.object$infectors)   #initialize nodes: will contain host carrying the node
    pbe1$v$nodetypes <- c(rep("s", sim.object$obs), rep("x", sim.object$Nsamples - sim.object$obs), 
                          rep("c", sim.object$Nsamples - 1), rep("t", sim.object$obs + 1))  #initialize nodes: will contain node type (sampling, additional sampling, coalescent)
  
    if(wh.bottleneck == "wide") {
      invisible(sapply(1:sim.object$obs, rewire_pullnodes_wide))
    } else {
      invisible(sapply(1:(sim.object$obs + 1), function(x){
        if (x == 1)
          pbe1$v$nodeparents[which(pbe1$v$nodehosts==0)] <- 0
        
        loosenodes <- which(pbe1$v$nodehosts == x & pbe1$v$nodeparents == -1)
        if (length(loosenodes) == 1){
          pbe1$v$nodeparents[loosenodes] <- 2*pbe1$d$nsamples+x-1
        } else {
          
          loosenodetimes <- pbe1$v$nodetimes[loosenodes]
          if (x == 1) {
            coalescenttimesnew <- sapply(loosenodetimes, function(t) runif(1, histtime, t))
            coalescenttimesnew <- coalescenttimesnew[-which(coalescenttimesnew == max(coalescenttimesnew))]
            coalescenttimesnew[which(coalescenttimesnew == min(coalescenttimesnew))] <- histtime
            
          } else {
            coalescenttimesnew <- sample_coaltimes(loosenodetimes, pbe1$v$inftimes[x], pbe1$p)
          }
          free_cnodes <- which(pbe1$v$nodeparents == -1 & pbe1$v$nodetypes == "c")
          
          coalnodes <- free_cnodes[1:(length(loosenodes)-1)]
          pbe1$v$nodehosts[coalnodes] <- x
          pbe1$v$nodeparents[coalnodes] <- 
            sapply(1:(length(loosenodes)-1), function(i){
              if (i == 1) return(2*pbe1$d$nsamples+x-1)
              else return(coalnodes[i-1])})
          pbe1$v$nodeparents[loosenodes] <- c(coalnodes,tail(coalnodes,1))
          pbe1$v$nodetimes[coalnodes] <- sort(coalescenttimesnew)
        }
      }))
    }  
    
  } else {
    pbe1$v$nodeparents <- rep(-1, 2 * sim.object$Nsamples + sim.object$obs - 1)  #initialize nodes: will contain parent node in phylotree
    pbe1$v$nodetimes <- c(sim.object$samtimes, sim.object$addsampletimes, 
                          rep(0, sim.object$Nsamples - 1), sim.object$inftimes)   #initialize nodes: will contain time of node
    pbe1$v$nodehosts <- c(1:sim.object$obs, sim.object$addsamplehosts, 
                          rep(-1, sim.object$Nsamples - 1), sim.object$infectors)   #initialize nodes: will contain host carrying the node
    pbe1$v$nodetypes <- c(rep("s", sim.object$obs), rep("x", sim.object$Nsamples - sim.object$obs), 
                          rep("c", sim.object$Nsamples - 1), rep("t", sim.object$obs))  #initialize nodes: will contain node type (sampling, additional sampling, coalescent)
 
    if(wh.bottleneck == "wide") {
      invisible(sapply(1:sim.object$obs, rewire_pullnodes_wide))
    } else {
      invisible(sapply(0:sim.object$obs, rewire_pullnodes_complete))
    }
  }
  res <- as.list.environment(pbe1)$v
  return(res)
  
}

### simulate sequences given a phylogenetic tree
sim_sequences <- function (sim.object, mu, sequence.length) {
  with(sim.object,{
    ### simulate the mutations on the phylotree
    #number of mutations
    edgelengths <- nodetimes - c(0, nodetimes)[1 + nodeparents]
    edgelengths[edgelengths < 0] <- 0  #rounding errors close to 0
    nmutations <- rpois(1, mu * sequence.length * sum(edgelengths))
    #place mutations on edges, order by time of edge (end)
    mutedges <- sample(length(edgelengths), size = nmutations, replace = TRUE, prob = edgelengths)
    mutedges <- mutedges[order(nodetimes[mutedges])]
    #sample mutations: which locus, to which nucleotide
    mutsites <- sample(sequence.length, size = nmutations, replace = TRUE)
    mutsites <- match(mutsites, unique(mutsites))  #place mutations at front of sequence 
    mutnucl <- sample(4, size = nmutations, replace = TRUE)
    
    ### construct the strains from the simulation by going backwards
    ### through the phylotree from each tip and placing mutations
    nodestrains <- matrix(data = rep(sample(4, nmutations, replace = TRUE), each = Nsamples), nrow = Nsamples)
    for(i in 1:Nsamples) {
      currentedge <- i
      recentmutations <- rep(FALSE, nmutations) #keep more recent mutations on each locus
      while(nodeparents[currentedge] != 0) {
        nodestrains[i, mutsites[mutedges == currentedge & !recentmutations]] <-
          mutnucl[mutedges == currentedge & !recentmutations]
        recentmutations <- recentmutations | mutedges == currentedge
        currentedge <- nodeparents[currentedge]
      }
    }
    # place single unchanged acgt at front of sequence to force these at first positions in phyDat-object
    nodestrains <- cbind(matrix(data = rep(1:4, each = Nsamples), ncol = 4), nodestrains)
    
    nodestrains[nodestrains == 1] <- "a"
    nodestrains[nodestrains == 2] <- "c"
    nodestrains[nodestrains == 3] <- "g"
    nodestrains[nodestrains == 4] <- "t"
    
    rownames(nodestrains) <- 1:Nsamples
    
    # make phyDat-object and change attributes to get entire sequence, with single acgt at front removed
    nodestrains <- phangorn::as.phyDat(nodestrains)
    mutlocs <- sample(4, max(0, sequence.length - nmutations), replace = TRUE)
    attr(nodestrains, "index") <- c(attr(nodestrains, "index")[-(1:4)], mutlocs)
    attr(nodestrains, "weight")[1:4] <- attr(nodestrains, "weight")[1:4] + tabulate(mutlocs, 4) - 1
    
    return(
      within(sim.object,{
        sequences <- nodestrains
      })
    )
  }
  )
}

### merge simulated trees together in one list
merge_trees <- function(sim.object) {
  obs <- sum(unlist(lapply(sim.object, function(xx) return(xx$obs))))
  Nsamples <- sum(unlist(lapply(sim.object, function(xx) return(xx$Nsamples))))
  
  sample.data <- do.call(rbind, lapply(sim.object, function(xx){
    return(data.frame("samtimes"=xx$samtimes,
                      "inftimes"=xx$inftimes,
                      "old"=1:length(xx$inftimes)))
  }))
  sample.data$tree <- unlist(lapply(seq_along(sim.object), function(i)
    return(rep(i, length(sim.object[[i]]$inftimes)))))
  merged <- order(sample.data$samtimes)
  sample.data <- sample.data[merged,]
  sample.data$merged <- merged
  sample.data$new <- 1:obs
  sample.data$infectors <- unlist(lapply(seq_along(sim.object), function(i) {
    inf <- sim.object[[i]]$infectors
    inf <- unlist(sapply(inf, function(xx){
      if (xx==0) return(0)
      return(sample.data$new[sample.data$tree==i & sample.data$old == xx])
    }))
    return(inf)
  }))[sample.data$merged]
  
  addsample.data <- do.call(rbind,lapply(seq_along(sim.object), function(i){
    hosts <- sim.object[[i]]$addsamplehosts
    hosts <- unlist(sapply(hosts, function(xx){
      if (xx==0) return(0)
      return(sample.data$new[sample.data$tree==i & sample.data$old == xx])
    }))
    return(data.frame("times"=sim.object[[i]]$addsampletimes,
                      "hosts"=hosts))
  }))
  if (length(addsample.data$times)!=0){
    addsample.data <- addsample.data[order(addsample.data$hosts),]
  } else {
    addsample.data <- list("times" = numeric(0),
                           "hosts" = integer(0))
  }
  
  nodehosts <- unlist(lapply(seq_along(sim.object), function(i) {
    host <- sim.object[[i]]$nodehosts
    host <- unlist(sapply(host, function(xx){
      if (xx==0) return(0)
      return(sample.data$new[sample.data$tree==i & sample.data$old == xx])
    }))
    return(host)
  }))
  
  node.data <- do.call(rbind, lapply(sim.object, function(xx) 
    return(data.frame("times"=xx$nodetimes, 
                      "type"=factor(xx$nodetypes,levels=c("s","x","c","t")),
                      "old"=1:length(xx$nodetypes)))))
  tree <- c()
  for (i in seq_along(sim.object)){
    tree <- c(tree, rep(i, length(sim.object[[i]]$nodetypes)))
  }
  node.data$tree <- tree
  node.data$hosts <- nodehosts
  node.order <- order(node.data$type, node.data$hosts, node.data$times)
  node.data <- node.data[node.order,]
  node.data$merged <- node.order
  node.data$new <- 1:nrow(node.data)
  node.data$parent <- unlist(lapply(seq_along(sim.object), function(i){
    xx <- sim.object[[i]]
    parent <- c()
    for (j in 1:length(xx$nodeparents)){
      newparent <- node.data$new[node.data$tree==i & node.data$old == xx$nodeparents[j]]
      if (length(newparent)==0) parent <- c(parent, 0)
      else parent <- c(parent,newparent)
    }
    return(parent)
  }))[node.order]
  
  sequences <- do.call(cbind,lapply(sim.object, function(xx) return(as.data.frame(xx$sequences))))
  seq.data <- do.call(rbind,lapply(sim.object, function(xx) 
    return(data.frame("times"=xx$nodetimes,
                      "type"=xx$nodetypes))))
  seq.data <- seq.data[seq.data$type == "s" | seq.data$type == "x",]
  seq.order <- order(seq.data$type,seq.data$times)
  sequences <- sequences[,seq.order,drop=FALSE]
  colnames(sequences) <- 1:ncol(sequences)
  if (Nsamples == 1){
    sequences <- cbind(sequences, sequences)
    colnames(sequences) <- 1:2
    
    sequences <- as.phyDat(sequences)
    sequences <- subset(sequences, subset = 1)
  } else {
    sequences <- as.phyDat(sequences)
  }
  
  return(list("obs"=obs,
              "infectors"=sample.data$infectors,
              "inftimes"=sample.data$inftimes,
              "samtimes"=sample.data$samtimes,
              "addsampletimes"=addsample.data$times,
              "addsamplehosts"=addsample.data$hosts,
              "Nsamples"=Nsamples,
              "nodeparents"=node.data$parent,
              "nodetimes"=node.data$times,
              "nodehosts"=node.data$hosts,
              "nodetypes"=as.character(node.data$type),
              "sequences"=sequences))
}


nthsample <- function(sim.object) {
  with(sim.object, {
    nth <- rep(0, Nsamples)
    sapply(1:obs, function(x) suppressWarnings(nth[which(nodehosts[nodetypes %in% c("s", "x")] == x)] <<- 0:Nsamples))
    return(nth)
  })
}
