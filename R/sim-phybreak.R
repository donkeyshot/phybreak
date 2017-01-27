### simulation of outbreaks in obkData format, according to the phybreak-model ###



#' Outbreak simulation.
#' 
#' Simulate outbreaks of class \code{'phybreakdata'}, with the outbreak model of \pkg{phybreak}.
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
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' simulation <- sim.phybreak()
#' @export
sim.phybreak <- function(obsize = 50, popsize = NA, Ngenes = 1,
                         R0 = 1.5, shape.gen = 10, mean.gen = 1, 
                         shape.sample = 10, mean.sample = 1,
                         wh.model = 3, wh.slope = 1, 
                         mu = 0.0001, sequence.length = rep(10000,Ngenes), output.class = c("list", "obkData"),
                         simmap = FALSE, reassortment = FALSE, reassortment.prob = reassortment*0.1) {
    
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
  if ( length(sequence.length) != Ngenes ) {
    stop("ngenes not equal to the given number of sequence lengths")
  }
  if (!is.logical(reassortment)) {
    stop("reassortment should be TRUE or FALSE")
  }
  if (reassortment.prob < 0 | reassortment.prob > 1){
    stop("reassortment probabilty should be between 0 and 1")
  }

  ### simulate step by step
  if(is.na(obsize)) {
    res <- .sim.outbreak(popsize, R0, shape.gen, mean.gen, shape.sample, mean.sample)
    obsize <- res$obs
    if(obsize == 1) stop(c("Outbreak size = 1"))
  } else {
    # Simulate outbreak for an observed size of cases, obsize.
    res <- .sim.outbreak.size(obsize, popsize, R0, shape.gen, mean.gen,
                              shape.sample, mean.sample) 
  }
  
  if(reassortment == TRUE & Ngenes>1){
    reassortment <- sample(c(0,1), obsize, replace = TRUE, prob = c(1-reassortment.prob, reassortment.prob))
  } else {
    reassortment <- rep(0, obsize)
  }
  
  # Simulate phylotree given a transmission tree
  res <- .sim.phylotree(res, wh.model, wh.slope, Ngenes, reassortment)  
  res <- .sim.sequences(res, mu, sequence.length, Ngenes) 
  hostnames <- paste0("host.", 1:obsize)
  
  ### make an obkData object 
  treesout <- vector('list',Ngenes)
  for(i in 1:Ngenes) {treesout[[i]] <-phybreak2phylo(vars = res, samplenames = hostnames, simmap = simmap, gene = i)}
  class(treesout) <- "multiPhylo" 
  names(treesout) <-  paste("Gene",1:Ngenes, sep = "" )
  
  sampletimes <- res$sampletimes
  names(sampletimes) <- hostnames
  infectiontimes <- res$infectiontimes
  names(infectiontimes) <- hostnames
  infectors <- c("index", hostnames)[1 + res$infector]
  names(infectors) <- hostnames
  seqs <- res$SNPlist
  names(seqs) <-  paste("SNPs_Gene",1:Ngenes, sep = "" )
  
  if(output.class == "obkData") {
    treesout <- converttreesout(treesout,Ngenes)
    toreturn <- new("obkData",
                    individuals = data.frame(
                      infector = infectors,
                      date = as.Date(infectiontimes, origin = "2000-01-01"),
                      row.names = hostnames),
                    dna = list(SNPs = seqs[[1]] ), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                    dna.individualID = hostnames, trees = treesout,
                    reassortment = reassortment)
    toreturn@dna@dna <- seqs
    toreturn@trees <- treesout
  } else {
    toreturn <- list(
      sequences = seqs,
      sample.times = sampletimes,
      sim.infection.times = infectiontimes,
      sim.infectors = infectors,
      sim.tree = treesout,
      reassortment = reassortment
    )
    class(toreturn) <- "phybreakdata"
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
    # when does currentID make infectious contacts?
    # reverse sorting so that with double contacts, the earliest will be used last
    whencontacts <- sort(inftimes[currentID] + rgamma(nrcontacts[currentID], aG, aG/mG),decreasing = TRUE)
    
    # who are these contacts made with?
    whocontacted <- sample(Npop, nrcontacts[currentID], replace = TRUE)
    
    # are these contacts successful, i.e. earlier than the existing contacts with these hosts?
    successful <- whencontacts < inftimes[whocontacted]
    
    # change infectors and infection times of successful contactees
    sources[whocontacted[successful]] <- currentID
    inftimes[whocontacted[successful]] <- whencontacts[successful]
    
    # go to next infected host in line
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
.sim.phylotree <- function (sim.object, wh.model, wh.slope, Ngenes, reassortment) {
  with(
    sim.object,
    {
      nodeparents <- matrix(rep( rep(0,3*obs - 1), Ngenes), nrow = Ngenes )  # initialize nodes: will containsparent node in phylotree
      nodetimes <- nodeparents        # initialize nodes: will contain time of node
      nodehosts <- nodeparents[1, ]   # initialize nodes: will contain host carrying the node
      nodetypes <- c(rep("s",obs), rep("c",obs - 1),
                     rep("t",obs))    # initialize nodes: will contain node type (sampling, coalescent, transmission)
      
      nodetimes[, 1:obs] <- matrix(rep(sampletimes,Ngenes), nrow = Ngenes, byrow = TRUE)                 # sampling nodes at sample times
      nodetimes[, 1:obs + 2*obs - 1] <- matrix(rep(infectiontimes,Ngenes), nrow = Ngenes, byrow = TRUE)  # transmission nodes at infection times
      
      nextcoalnode <- obs + 1  # needed in next loop: which is the next available coalescence node
      
      ## distribute nodes over the hosts
      nodehosts[1:obs] <- 1:obs
      nodehosts[(obs+1):(2*obs-1)] <- sort(infectors)[-1]
      nodehosts[(2*obs):(3*obs-1)] <- infectors
      
      for(i in 1:obs) {
        iees <- which(i == infectors)  # infectees of host i (transmission nodes in host i)
        nodehosts[i] <- i              # assign sample node
        while(length(iees) > 0) {      # per infectee...
          nodehosts[iees[1] + 2*obs - 1] <- i  # ...assign transmission node
          nodehosts[nextcoalnode] <- i         # ...assign coalescence node
          nextcoalnode <- nextcoalnode + 1     # next coalescence node
          iees <- iees[-1]                     # infectee done with
        }
      }
      
      ## sample the times of the coalescence nodes
      for(i in 1:obs) {
        if (reassortment[i] == 1){
          for(j in 1:Ngenes) {
            nodetimes[j, nodehosts == i & nodetypes == "c"] <-   # change the times of the coalescence nodes in host i...
              nodetimes[j, i + 2*obs - 1] +                      # ...to the infection time +
              .samplecoaltimes(nodetimes[j, nodehosts == i & nodetypes != "c"] - nodetimes[j, i + 2*obs - 1],
                               wh.model, wh.slope)               # ...sampled coalescence times
            
            nodeparents[j, nodehosts == i] <-     #change the parent nodes of all nodes in host i...
              .sampletopology(which(nodehosts == i), nodetimes[j, nodehosts == i], nodetypes[nodehosts == i], i + 2*obs - 1, wh.model)
            #...to a correct topology, randomized where possible
          } 
        } else { 
          coalt <-  nodetimes[1, i + 2*obs - 1] + 
            .samplecoaltimes(nodetimes[1, nodehosts == i & nodetypes != "c"] - nodetimes[1, i + 2*obs - 1],
                             wh.model, wh.slope)  
          nodetimes[, nodehosts == i & nodetypes == "c"] <- matrix(rep(coalt,Ngenes),nrow = Ngenes, byrow = TRUE)
          
          nodep <-  matrix(rep(.sampletopology(which(nodehosts == i), nodetimes[1, nodehosts == i],
                                               nodetypes[nodehosts == i], i + 2*obs - 1, wh.model), Ngenes),
                           byrow = TRUE, nrow = Ngenes)
          nodeparents[, nodehosts == i] <- nodep
        }
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
.sim.sequences <- function (sim.object, mu, sequence.length, Ngenes) {
  if (Ngenes > 1) SNPList <- list() else SNPList <- c()
  with(sim.object,{
    for (j in 1:Ngenes) {
      ### simulate the mutations on the phylotree
      # number of mutations
      edgelengths <- nodetimes[j, ] - c(0,nodetimes[j, ])[1+nodeparents[j, ]]
      edgelengths[edgelengths < 0] <- 0  # rounding errors close to 0
      nmutations <- rpois(1, mu * sequence.length[j] * sum(edgelengths))
      if(nmutations == 0) {nmutations <- 1}
      # place mutations on edges, order by time of edge (end)
      mutedges <- sample(3*obs-1, size = nmutations, replace = TRUE, prob = edgelengths)
      mutedges <- mutedges[order(nodetimes[j, mutedges])]
      # sample mutations: which locus, to which nucleotide
      mutsites <- sample(sequence.length[j], size = nmutations, replace = TRUE)
      mutsites <- match(mutsites, unique(mutsites))  # place mutations at front of sequence
      mutnucl <- sample(4, size = nmutations, replace = TRUE)
      
      ### construct the strains from the simulation by going backwards
      ### through the phylotree from each tip and placing mutations
      nodestrains <- matrix(data = rep(sample(4, nmutations, replace = TRUE), each = obs), ncol = nmutations)
      for(i in 1:obs) {
        currentedge <- i                          
        recentmutations <- rep(FALSE, nmutations) # keep more recent mutations on each locus
        while(nodeparents[j,currentedge] != 0) {          
          nodestrains[i,mutsites[mutedges == currentedge & !recentmutations]] <- 
            mutnucl[mutedges == currentedge & !recentmutations]                  
          recentmutations <- recentmutations | mutedges == currentedge           
          currentedge <- nodeparents[j, currentedge]               
        }                                                                         
      }
      # place single unchanged acgt at front of sequence to force these at first positions in phyDat-object
      nodestrains <- cbind(matrix(data = rep(1:4, each = obs), ncol = 4), nodestrains)
      
      nodestrains[nodestrains == 1] <- "a"
      nodestrains[nodestrains == 2] <- "c"
      nodestrains[nodestrains == 3] <- "g"
      nodestrains[nodestrains == 4] <- "t"
      
      rownames(nodestrains) <- paste0("host.",1:obs)
      
      # make phyDat-object and change attributes to get entire sequence, with single acgt at front removed
      nodestrains <- phangorn::as.phyDat(nodestrains)
      mutlocs <- sample(4, sequence.length[j] - nmutations, replace = TRUE)
      attr(nodestrains, "index") <- c(attr(nodestrains, "index")[-(1:4)], mutlocs)
      attr(nodestrains, "weight")[1:4] <- attr(nodestrains, "weight")[1:4] + tabulate(mutlocs, 4) - 1
      
      SNPList[[j]] <- ape::as.DNAbin(nodestrains) 
    }
    return(within(sim.object,{SNPlist <- SNPList}))
  })
}


### Tree information should be modified for multiple genes in order to build an obkData object.
### called by:
# sim.phybreak
converttreesout <- function(treesout, Ngenes){
  trees <- c()
  
  for (i in 1:Ngenes){
    trees$edge[[i]] <- treesout[[i]]$edge
    trees$edge.length[[i]] <- treesout[[i]]$edge.length
    trees$Nnode <- treesout[[1]]$Nnode
    trees$tip.label <- treesout[[1]]$tip.label
    trees$root.edge[[i]] <- treesout[[i]]$root.edge
  }
  
  names(trees$edge) <-  paste("Gene",1:Ngenes, sep = "" )
  names(trees$edge.length) <-  paste("Gene",1:Ngenes, sep = "" )
  names(trees$root.edge) <-  paste("Gene",1:Ngenes, sep = "" )
  class(trees) <- "multiPhylo"
  
  return(trees)
}
### Extract a tree for gene x from either an Ngenes>1 obkData or phybreak object
extractTree <- function(data.object, gene){
  
  if (class(data.object) == "phybreakdata") {
    tree <- data.object
    tree$d$sequences <- data.object$d$sequences[[gene]]
    tree$v$nodetimes <- data.object$v$nodetimes[gene, ]
    tree$v$nodeparents <- data.object$v$nodeparents[gene,]
    tree$s$nodetimes <- data.object$s$nodetimes[[gene]]
    tree$s$nodeparents <- data.object$s$nodeparents[[gene]]
    tree$s$logLikseq <- data.object$s$logLikseq[genes, ]
    tree$d$nSNPs <- data.object$d$nSNPs[gene]
  } else if (class(data.object) == "obkData"){
    tree <- data.object
    edge <- data.object@trees$edge[[gene]]
    edge.length <-data.object@trees$edge.length[[gene]]
    root.edge <-  data.object@trees$root.edge[[gene]]
    DNA <- data.object@dna@dna[[gene]]
    
    # Build tree for a selected gene
    trees <- list()
    trees$edge <- edge
    trees$edge.length <- edge.length  
    trees$Nnode <- data.object@trees$Nnode
    trees$tip.label <- data.object@trees$tip.label
    trees$root.edge <- root.edge
    class(trees) <- "phylo"
    class(trees) <- "multiPhylo"
    
    tree@trees <- trees
    tree@dna@dna <- NULL
    tree@dna@dna$SNPs <- DNA
  } else {
    stop("data object input is not of class 'obkData' or 'phybreakdata' ")
  }
  return(do)
}