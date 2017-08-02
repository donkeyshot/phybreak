### simulation of outbreaks in obkData format, according to the phybreak-model ###



#' Outbreak simulation.
#' 
#' Simulate outbreaks of class \code{'phybreakdata'}, with the outbreak model of \pkg{phybreak}.
#' 
#' @param obsize The outbreak size (number of cases) to obtain. If \code{obsize = NA}, \code{popsize} should be provided.
#' @param popsize The population size in which to simulate. If it is not defined (default), 
#'   an optimal population size will be chosen based on R0 and obsize. Be aware that choosing a \code{popsize} and
#'   an \code{obsize} can severely increase the simulation time, depending on \code{R0}.
#' @param samplesperhost Number of samples to be taken per host, either a vector or a single number.
#' @param R0 The basic reproduction ratio used for simulation. The offspring distribution is Poisson.
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
#'     \item "single": effective size = 0, so coalescence occurs 'just before' transmission in the infector (strict bottleneck)
#'     \item "infinite": effective size = Inf, with strict bottleneck, so coalescence occurs 'just after' transmission in the infectee
#'     \item "linear": effective size at time t after infection = \code{wh.slope * t} (strict bottleneck)
#'     \item "exponential": effective size at time t after infection = \code{wh.level * exp(wh.exponent * t)} (loose bottleneck)
#'     \item "constant": effective size = wh.level (loose bottleneck)
#'   }
#' @param wh.slope Initial value for the within-host slope, used if \code{wh.model = "linear"}.
#' @param wh.exponent Initial value for the within-host exponent, used if \code{wh.model = "exponential"}
#' @param wh.level Initial value for the within-host effective pathogen size at transmission, used if 
#'   \code{wh.model = "exponential"} or if \code{wh.model = "constant"}
#' @param mu Expected number of mutations per nucleotide per unit of time along each lineage. 
#' @param sequence.length Number of available nucleotides for mutations.
#' @param output.class Class of the simulation output. If package \pkg{OutbreakTools} is available, it is possible to choose
#'  class \code{'obkData'}
#' @param ... If arguments from previous versions of this function are used, they may be interpreted correctly through 
#'   this argument, but it is better to provide the correct argument names.
#' @return The simulation output, either as an object of class \code{'phybreakdata'} with sequences (class \code{'phyDat'}) and 
#'   sampling times (which would be the observations), and infection times, infectors, and phylogenetic tree 
#'   of class \code{\link[ape]{phylo}}; 
#'   or as an object of class \code{'obkData'} (package \pkg{OutbreakTools}), containing the outbreak data in the following slots:
#'   \describe{
#'     \item{individuals}{a \code{data.frame} with individual labels as row names, a vector \code{infector},
#'       and a vector \code{date} containing the infection times (starting 01-01-2000)
#'     }
#'     \item{dna}{an object of class \code{'obkSequences'}, with SNP data in \code{dna} and sampling times
#'       in \code{meta$date}
#'     }
#'     \item{trees}{an object of class \code{\link[ape]{multiphylo}}, containing a single tree of class \code{\link[ape]{phylo}}}
#'   }
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' simulation <- sim.phybreak()
#' @export
sim.phybreak <- function(obsize = 50, popsize = NA, samplesperhost = 1,
                         R0 = 1.5, gen.shape = 10, gen.mean = 1, 
                         sample.shape = 10, sample.mean = 1,
                         additionalsampledelay = 0,
                         wh.model = "linear", wh.slope = 1, wh.exponent = 1, wh.level = 0.1,
                         mu = 0.0001, sequence.length = 10000, output.class = c("phybreakdata", "obkData"), ...) {
  ### parameter name compatibility 
  old_arguments <- list(...)
  if(exists("oldarguments$shape.gen")) gen.shape <- oldarguments$shape.gen
  if(exists("oldarguments$mean.gen")) gen.mean <- oldarguments$mean.gen
  if(exists("oldarguments$shape.sample")) sample.shape <- oldarguments$shape.sample
  if(exists("oldarguments$mean.sample")) sample.mean <- oldarguments$sample.mean

  
  ### tests
  output.class <- match.arg(output.class)
  if(output.class == "obkData" && !("OutbreakTools" %in% .packages(TRUE))) {
    warning("package 'OutbreakTools' not installed: output.class is \"phybreakdata\"")
    output.class <- "phybreakdata"
  }
  if(output.class == "obkData" && !("OutbreakTools" %in% .packages(FALSE))) {
#    warning("package 'OutbreakTools' is not attached")
    requireNamespace("OutbreakTools")
  }
  if(all(is.na(c(obsize, popsize)))) {
    stop("give an outbreak size (obsize) and/or a population size (popsize)")
  }
  if(all(!is.na(c(obsize, popsize)))) {
    warning("giving both an outbreak size (obsize) and a population size (popsize) can take a very long simulation time",
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
  
  ### simulate step by step
  if(is.na(obsize)) {
     res <- .sim.outbreak(popsize, R0, gen.shape, gen.mean,
                                      sample.shape, sample.mean)
     obsize <- res$obs
     if(obsize == 1) return(c("Outbreak size = 1"))
  } else {
    res <- .sim.outbreak.size(obsize, popsize, R0, gen.shape, gen.mean,
                              sample.shape, sample.mean)
  }
  if(any(samplesperhost < 1)) stop("samplesperhost should be positive")
  if(any(additionalsampledelay < 0)) stop("additionalsampledelay cannot be negative")
  res <- .sim.additionalsamples(res, samplesperhost, additionalsampledelay)
  res <- .sim.phylotree(res, wh.model, wh.slope, wh.exponent, wh.level, sample.mean)
  res <- .sim.sequences(res, mu, sequence.length)
  hostnames <- paste0("host.", 1:obsize)
  samplenames <- paste0("sample.", res$nodehosts[1:res$Nsamples], ".", nthsample(res))
  names(res$sequences) <- samplenames
  
  ### make a phylo tree
  treesout <- vector('list',1)
  treesout[[1]] <- phybreak2phylo(vars = res, samplenames = samplenames, simmap = FALSE)
  class(treesout) <- "multiPhylo"
  

  if(output.class == "obkData") {
    toreturn <- with(res, new("obkData",
                    individuals = data.frame(
                      infector = c("index", hostnames)[1 + infectors],
                      date = as.Date(inftimes, origin = "2000-01-01"),
                      row.names = hostnames),
                    dna = list(SNPs = ape::as.DNAbin(sequences)), 
                    dna.individualID = hostnames[c(1:obs, addsamplehosts)], 
                    dna.date = as.Date(c(samtimes, addsampletimes), origin = "2000-01-01"),
                    sample = samplenames,
                    trees = treesout))
  } else {
    toreturn <- with(res,
                     phybreakdata(
      sequences = sequences,
      sample.times = c(samtimes, addsampletimes),
      sample.names = samplenames,
      host.names = hostnames[c(1:obs, addsamplehosts)],
      sim.infection.times = inftimes,
      sim.infectors = infectors,
      sim.tree = treesout[[1]]
    ))
  }
  return(toreturn)
}


### simulate an outbreak of a particular size by repeating
### simulations until obsize is obtained
### called by:
# sim.phybreak
### calls:
# .sim.outbreak
.sim.outbreak.size <- function(obsize, Npop, R0, aG, mG, aS, mS) {
  if(is.na(Npop)) {
    Npop <- obsize
    while(1 - obsize/Npop < exp(-R0* obsize/Npop)) {Npop <- Npop + 1}
  } 
  
  sim <- .sim.outbreak(Npop, R0, aG, mG, aS, mS)
  
  while(sim$obs != obsize) {
    sim <- .sim.outbreak(Npop, R0, aG, mG, aS, mS)
  }
  
  return(sim)
}

### simulate an outbreak
### called by:
# sim.phybreak
# .sim.outbreak.size
.sim.outbreak <- function(Npop, R0, aG, mG, aS, mS) {
  ### initialize
  inftimes <- c(0, rep(10000, Npop-1))
  sources <- rep(0,Npop)
  nrcontacts <- rpois(Npop, R0)
  nth.infection <- 1
  currentID <- 1
  
  ### by order of infection, sample secondary infections
  # currentID is the infected host under consideration
  while(nth.infection <= Npop & inftimes[currentID] != 10000) {
    #when does currentID make infectious contacts?
    #reverse sorting so that with double contacts, the earliest will be used last
    whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
    
    #who are these contacts made with?
    whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE)
    
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
  samtimes <- inftimes + rgamma(Npop, aS, aS/mS)
  
  ### order hosts by sampling times, and renumber hostIDs
  ### so that the uninfected hosts can be discarded
  orderbysamtimes <- order(samtimes)
  sources <- sources[orderbysamtimes]
  infectors <- match(sources,orderbysamtimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamtimes]
  samtimes <- samtimes[orderbysamtimes]
  
  ### return the outbreak
  return(
    list(
      obs = obs,
      samtimes = samtimes[1:obs],
      inftimes = inftimes[1:obs],
      infectors = infectors
    )
  )
}


### simulate additional samples in a transmission tree
.sim.additionalsamples <- function(sim.object, samperh, addsamdelay) {
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
### called by:
# sim.phybreak
### calls:
# .samplecoaltimes
# .sampletopology
.sim.phylotree <- function (sim.object, wh.model, wh.slope, wh.exponent, wh.level, sample.mean) {
  simenv <- new.env()
  list2env(list(v = sim.object, 
                p = list(wh.model = wh.model, wh.slope = wh.slope, wh.exponent = wh.exponent,
                         wh.level = wh.level, sample.mean = sample.mean)), simenv)
  simenv$v$nodeparents <- c(0, (sim.object$Nsamples + 1):(2*sim.object$Nsamples - 1), 
                            rep(-1, sim.object$Nsamples - 1))  #initialize nodes: will contain parent node in phylotree
  simenv$v$nodetimes <- c(sim.object$samtimes, sim.object$addsampletimes, sim.object$samtimes[-1], 
                         sim.object$addsampletimes)   #initialize nodes: will contain time of node
  simenv$v$nodehosts <- c(1:sim.object$obs, sim.object$addsamplehosts, 
                         2:sim.object$obs, sim.object$addsamplehosts)   #initialize nodes: will contain host carrying the node
  simenv$v$nodetypes <- c(rep("s", sim.object$obs), rep("x", sim.object$Nsamples - sim.object$obs), 
                         rep("c", sim.object$Nsamples - 1))  #initialize nodes: will contain node type (sampling, additional sampling, coalescent)
  invisible(sapply(1:sim.object$obs, rewire_buildminitree, phybreakenv = simenv))
  res <- as.list.environment(simenv)$v
  return(res)
  
}


### simulate sequences given a phylogenetic tree
### called by:
# sim.phybreak
.sim.sequences <- function (sim.object, mu, sequence.length) {
  with(sim.object,{
    ### simulate the mutations on the phylotree
    #number of mutations
    edgelengths <- nodetimes - c(0, nodetimes)[1 + nodeparents]
    edgelengths[edgelengths < 0] <- 0  #rounding errors close to 0
    nmutations <- rpois(1, mu * sequence.length * sum(edgelengths))
    #place mutations on edges, order by time of edge (end)
    mutedges <- sample(2 * Nsamples - 1, size = nmutations, replace = TRUE, prob = edgelengths)
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


nthsample <- function(sim.object) {
  with(sim.object, {
    nth <- rep(0, Nsamples)
    sapply(1:obs, function(x) suppressWarnings(nth[which(nodehosts[nodetypes %in% c("s", "x")] == x)] <<- 0:Nsamples))
    return(nth)
  })
}
