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
.sim.phylotree <- function (sim.object, wh.model, slope) {
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
                           wh.model, slope)  #...sampled coalescence times
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
.sim.sequences <- function (sim.object, mutrate, sequence.length) {
  with(sim.object,{
    ### simulate the mutations on the phylotree
    #number of mutations
    edgelengths <- nodetimes - c(0,nodetimes)[1+nodeparents]
    edgelengths[edgelengths < 0] <- 0  #rounding errors close to 0
    nmutations <- rpois(1, mutrate * sequence.length * sum(edgelengths))
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
        SNPlist <- as.DNAbin(nodestrains)
      })
    )
  }
  )
}
