### simulation of outbreaks in obkData format, according to the phybreak-model ###

### calls:
# .sim.outbreak.size
# .sim.outbreak
# .sim.phylotree
# .sim.sequences
# .makephylo.phybreak
### required packages:
# OutbreakTools
# phangorn


#' Outbreak simulation.
#' 
#' Simulate outbreaks of class \code{'obkData'} (package \pkg{OutbreakTools}), with the outbreak model of \pkg{phybreak}.
#' 
#' @param obsize The outbreak size (number of cases) to obtain. If \code{obsize = NA}, \code{popsize} should be provided.
#' @param popsize The population size in which to simulate. If it is not defined (default), 
#'   an optimal population size will be chosen based on R0 and obsize. Be aware that choosing a \code{popsize} and
#'   an \code{obsize} can severely increase the simulation time, depending on \code{R0}.
#' @param R0 The basic reproduction ratio used for simulation. The offspring distribution is Poisson.
#' @param shape.gen The shape parameter of the gamma-distributed generation interval.
#' @param mean.gen The mean generation interval.
#' @param shape.sample The shape parameter of the gamma-distributed sampling interval.
#' @param mean.sample The mean sampling interval.
#' @param wh.model The model for within-host pathogen dynamics (effective pathogen population size = 
#'   N*gE = actual population size * pathogen generation time), used to simulate coalescence events. Options are:
#'   \enumerate{
#'     \item Effective size = 0, so coalescence occurs 'just before' transmission, in the infector
#'     \item Effective size = Inf, so coalescence occurs 'just after' infection, in the infectee
#'     \item Effective size at time t after infection = wh.slope * t
#'   }
#' @param wh.slope Within-host increase of effective population size, used if \code{wh.model = 3}.
#' @param mu Expected number of mutations per nucleotide per unit of time along each lineage. 
#' @param sequence.length Number of available nucleotides for mutations.
#' @param output.class Class of the simulation output. If package \pkg{OutbreakTools} is available, it is possible to choose
#'  class \code{'obkData'}
#' @return The simulation output, either in a list with sequences (class \code{'phyDat'}) and sampling times (which would be
#'   the observations), and infection times and infectors; or as an object of class \code{'obkData'} (package \pkg{OutbreakTools}), 
#'   containing the outbreak data in the following slots:
#'   \describe{
#'     \item{individuals}{a \code{data.frame} with individual labels as row names, a vector \code{infector},
#'       and a vector \code{date} containing the infection times (starting 01-01-2000)
#'     }
#'     \item{dna}{an object of class \code{'obkSequences'}, with SNP data in \code{dna} and sampling times
#'       in \code{meta$date}
#'     }
#'     \item{trees}{an object of class \code{\link[ape]{multiphylo}}, containing a single tree of class \code{phylo}}
#'   }
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' simulation <- sim.phybreak()
#' @export
sim.phybreak <- function(obsize = 50, popsize = NA, 
                         R0 = 1.5, shape.gen = 10, mean.gen = 1, 
                         shape.sample = 10, mean.sample = 1,
                         wh.model = 3, wh.slope = 1, 
                         mu = 0.0001, sequence.length = 10000, output.class = c("list", "obkData")) {
  ### tests
  output.class <- output.class[output.class %in% c("list", "obkData")][1]
  if(is.na(output.class)) {
    stop("output.class should be \"list\" or \"obkData\"")
  }
  if(output.class == "obkData" & !("OutbreakTools" %in% .packages(TRUE))) {
    warning("package 'OutbreakTools' not installed: output.class is \"list\"")
    output.class <- "list"
  }
  if(output.class == "obkData" & !("OutbreakTools" %in% .packages(FALSE))) {
#    warning("package 'OutbreakTools' is not attached")
    requireNamespace("OutbreakTools")
  }
  if(all(is.na(c(obsize,popsize)))) {
    stop("give an outbreak size (obsize) and/or a population size (popsize)")
  }
  if(all(!is.na(c(obsize,popsize)))) {
    warning("giving both an outbreak size (obsize) and a population size (popsize) can take a very long simulation time",
            immediate. = TRUE)
  }
  if(all(!is.na(c(obsize,popsize))) && obsize > popsize) {
    stop("outbreak size (obsize) cannot be larger than population size (popsize)")
  }
  if(R0 <= 1) {
    stop("R0 should be larger than 1")
  }
  if(any(c(shape.gen, mean.gen, shape.sample, mean.sample, wh.slope, mu) <= 0)) {
    stop("parameter values should be positive")
  }
  
  ### simulate step by step
  if(is.na(obsize)) {
     res <- .sim.outbreak(popsize, R0, shape.gen, mean.gen,
                                      shape.sample, mean.sample)
     obsize <- res$obs
     if(obsize == 1) return(c("Outbreak size = 1"))
  } else {
    res <- .sim.outbreak.size(obsize, popsize, R0, shape.gen, mean.gen,
                              shape.sample, mean.sample)
  }
  res <- .sim.phylotree(res, wh.model, wh.slope)
  res <- .sim.sequences(res, mu, sequence.length)

  ### make an obkData object
  treesout <- vector('list',1)
  treesout[[1]] <- .makephylo.phybreak(res$nodetimes, res$nodeparents, 1:obsize)
  class(treesout) <- "multiPhylo"
  
  if(output.class == "obkData") {
    toreturn <- new("obkData",
                    individuals = data.frame(
                      infector = res$infectors,
                      date = as.Date(res$infectiontimes, origin = "2000-01-01"),
                      row.names = 1:res$obs),
                    dna = list(SNPs = res$SNPlist), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                    dna.individualID = 1:res$obs, trees = treesout)
  } else {
    sampletimes <- res$sampletimes
    names(sampletimes) <- paste0("host.", 1:res$obs)
    seqs <- phangorn::as.phyDat(res$SNPlist)
    names(seqs) <- paste0("host.", 1:res$obs)
    infectiontimes <- res$infectiontimes
    names(infectiontimes) <- paste0("host.", 1:res$obs)
    infectors <- c("index", paste0("host.", 1:res$obs))[1 + res$infector]
    toreturn <- list(
      sequences = seqs,
      sample.times = sampletimes,
      infection.times = infectiontimes,
      infectors = infectors
    )
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
  samptimes <- inftimes + rgamma(Npop, aS, aS/mS)
  
  ### order hosts by sampling times, and renumber hostIDs
  ### so that the uninfected hosts can be discarded
  orderbysamptimes <- order(samptimes)
  sources <- sources[orderbysamptimes]
  infectors <- match(sources,orderbysamptimes)[1:obs]
  infectors[is.na(infectors)] <- 0
  inftimes <- inftimes[orderbysamptimes]
  samptimes <- samptimes[orderbysamptimes]
  
  ### return the outbreak
  return(
    list(
      obs = obs,
      sampletimes = samptimes[1:obs],
      infectiontimes = inftimes[1:obs],
      infectors = infectors
    )
  )
}


### simulate a phylogenetic tree given a transmission tree
### called by:
# sim.phybreak
### calls:
# .samplecoaltimes
# .sampletopology
.sim.phylotree <- function (sim.object, wh.model, wh.slope) {
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
                           wh.model, wh.slope)  #...sampled coalescence times
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


### simulate sequences given a phylogenetic tree
### called by:
# sim.phybreak
.sim.sequences <- function (sim.object, mu, sequence.length) {
  with(sim.object,{
    ### simulate the mutations on the phylotree
    #number of mutations
    edgelengths <- nodetimes - c(0,nodetimes)[1+nodeparents]
    edgelengths[edgelengths < 0] <- 0  #rounding errors close to 0
    nmutations <- rpois(1, mu * sequence.length * sum(edgelengths))
    #place mutations on edges, order by time of edge (end)
    mutedges <- sample(3*obs-1, size = nmutations, replace = TRUE, prob = edgelengths)
    mutedges <- mutedges[order(nodetimes[mutedges])]
    #sample mutations: which locus, to which nucleotide
    mutsites <- sample(sequence.length, size = nmutations, replace = TRUE)
    mutsites <- match(mutsites, unique(mutsites))  #place mutations at front of sequence 
    mutnucl <- sample(4, size = nmutations, replace = TRUE)
    
    ### construct the strains from the simulation by going backwards
    ### through the phylotree from each tip and placing mutations
    nodestrains <- matrix(data = rep(1, sequence.length * obs), ncol = sequence.length)
    for(i in 1:obs) {
      currentedge <- i
      recentmutations <- rep(FALSE, nmutations) #keep more recent mutations on each locus
      while(nodeparents[currentedge] != 0) {
        nodestrains[i,mutsites[mutedges == currentedge & !recentmutations]] <-
          mutnucl[mutedges == currentedge & !recentmutations]
        recentmutations <- recentmutations | mutedges == currentedge
        currentedge <- nodeparents[currentedge]
      }
      
    }
    nodestrains[nodestrains == 1] <- "a"
    nodestrains[nodestrains == 2] <- "c"
    nodestrains[nodestrains == 3] <- "g"
    nodestrains[nodestrains == 4] <- "t"
    
    rownames(nodestrains) <- 1:obs
    
    return(
      within(sim.object,{
        SNPlist <- ape::as.DNAbin(nodestrains)
      })
    )
  }
  )
}
