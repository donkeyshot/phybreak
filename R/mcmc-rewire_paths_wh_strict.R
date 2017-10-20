rewire_pathA_wh_strict <- function() {
  ### First, dismantle minitree
  # edges entering hostID, with endtimes
  edgesin <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c")
  edgeintimes <- pbe1$v$nodetimes[edgesin]

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # coalescent nodes before and within hostID
  coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")

  # dismantle topology, move transmission node
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop

  ### Second, change transmission tree
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop

  ### Third, rebuild minitree
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(edgeintimes, pbe1$tinf.prop, pbe1$p)

  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]

  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)

  # change minitree in hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes
}

rewire_pathB_wh_strict <- function() {
  ### Identify new index
  newindex <- which(pbe1$v$inftimes == sort(pbe1$v$inftimes)[2])

  ### First, dismantle minitree
  # transmission node of new index
  transnode_ni <- 2*pbe1$d$nsamples - 1 + newindex

  # edges entering hostID, with endtimes (excluding t from new index)
  edgesin <- setdiff(which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"),
                     transnode_ni)
  edgeintimes <- pbe1$v$nodetimes[edgesin]

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # coalescent nodes in hostID
  coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")

  # dismantle topology, move transmission and bottleneck nodes
  pbe1$v$nodehosts[coalnodes] <- -1L
  pbe1$v$nodehosts[transnode_ni] <- 0L
  pbe1$v$nodehosts[transnode] <- pbe1$infector.proposed.ID
  pbe1$v$nodeparents[c(edgesin, coalnodes, transnode, transnode_ni)] <- -1
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop

  ### Second, change transmission tree
  pbe1$v$infectors[pbe1$hostID] <- pbe1$infector.proposed.ID
  pbe1$v$infectors[newindex] <- 0L
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop

  ### Third, rebuild minitree in hostID
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(edgeintimes, pbe1$tinf.prop, pbe1$p)
  coalnodes <- tail(coalnodes, -1)

  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]

  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)

  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- pbe1$hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes

  ### Fourth, add edges in infector and new index
  rewire_pullnodes_wh_strict(0)
  rewire_pullnodes_wh_strict(pbe1$infector.proposed.ID)
}

rewire_pathCF1_wh_strict <- function() {
  ### Identify new infector and old infector
  newinfector <- which(pbe1$v$infectors == pbe1$hostID)
  newinfector <- newinfector[which(pbe1$v$inftimes[newinfector] == min(pbe1$v$inftimes[newinfector]))]
  oldinfector <- pbe1$v$infectors[pbe1$hostID]

  ### First, dismantle minitree
  # edges entering new infector, with endtimes
  edgesin_ni <- which(pbe1$v$nodehosts == newinfector & pbe1$v$nodetypes != "c")
  edgeintimes_ni <- pbe1$v$nodetimes[edgesin_ni]

  # transmission node of new infector
  transnode_ni <- 2*pbe1$d$nsamples - 1 + newinfector

  # edges entering hostID, with endtimes, excluding t from new infector
  edgesin <- setdiff(which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c"),
                     transnode_ni)
  edgeintimes <- pbe1$v$nodetimes[edgesin]

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # all coalescent nodes in new infector and hostID
  coalnodes <- which((pbe1$v$nodehosts == newinfector | pbe1$v$nodehosts == pbe1$hostID) & pbe1$v$nodetypes == "c")

  # more coalescent nodes:
  if(oldinfector > 0) {
    # parent of transmission edge leaving hostID
    coalnodes <- c(take_cnode(transnode), coalnodes)
  }

  # dismantle topology, move transmission nodes
  pbe1$v$nodehosts[coalnodes] <- -1L
  pbe1$v$nodehosts[transnode_ni] <- oldinfector
  pbe1$v$nodehosts[transnode] <- newinfector
  pbe1$v$nodeparents[c(edgesin, edgesin_ni, coalnodes, transnode, transnode_ni)] <- -1L
  pbe1$v$nodetimes[c(transnode, transnode_ni)] <- pbe1$v$nodetimes[c(transnode_ni, transnode)]

  ### Second, change transmission tree
  pbe1$v$infectors[newinfector] <- oldinfector
  pbe1$v$infectors[pbe1$hostID] <- newinfector
  pbe1$v$inftimes[c(newinfector, pbe1$hostID)] <- pbe1$v$inftimes[c(pbe1$hostID, newinfector)]

  ### Third, rebuild minitree in hostID
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(edgeintimes, pbe1$v$inftimes[pbe1$hostID], pbe1$p)
  coalnodes_ni <- head(coalnodes, length(coalnodes) - length(newcoaltimes))
  coalnodes <- tail(coalnodes, length(newcoaltimes))

  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]

  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)

  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- pbe1$hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes

  ### Fourth, rebuild minitree in new infector
  # add edges from hostID to incoming edges
  edgesin_ni <- c(edgesin_ni, transnode)
  edgeintimes_ni <- pbe1$v$nodetimes[edgesin_ni]

  # times of coalescent events in new infector and bottleneck size, and distribute coalescent nodes over new infector and pre-newinfector
  newcoaltimes <- sample_coaltimes(edgeintimes_ni, pbe1$v$inftimes[newinfector], pbe1$p)
  coalnodes_ni <- tail(coalnodes_ni, length(newcoaltimes))

  # order all edges (transmission, sample, coalescent) by endtime in new infector
  nodeorder <- order(c(newcoaltimes, edgeintimes_ni))
  edgeend <- c(coalnodes_ni, edgesin_ni)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes_ni)[nodeorder]

  # sample topology of minitree within new infector
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes_ni)))[nodeorder],
                               transnode_ni)

  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- newinfector
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes

  ### Fifth, add edges before new infector
  rewire_pullnodes_wh_strict(oldinfector)
}

rewire_pathD_wh_strict <- function() {
  ### First, dismantle minitree
  # edges entering hostID, with endtimes
  edgesin <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c")
  edgeintimes <- pbe1$v$nodetimes[edgesin]

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # coalescent nodes of hostID
  coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")

  # more coalescent nodes:
  # parent of transmission edge leaving hostID
  coalnodes <- c(take_cnode(transnode), coalnodes)

  # new edges entering hostID, from old index
  newedgesin <- which(pbe1$v$nodehosts == 0)
  newedgeintimes <- pbe1$v$nodetimes[newedgesin]

  # dismantle topology, move transmission and bottleneck nodes
  pbe1$v$nodehosts[coalnodes] <- -1L
  pbe1$v$nodehosts[newedgesin] <- pbe1$hostID
  pbe1$v$nodehosts[transnode] <- 0L
  pbe1$v$nodeparents[c(edgesin, newedgesin, coalnodes, transnode)] <- -1L
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop

  ### Second, change transmission tree
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop
  pbe1$v$infectors[pbe1$v$infectors == 0] <- pbe1$hostID
  pbe1$v$infectors[pbe1$hostID] <- 0L

  ### Third, rebuild minitree in hostID
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(c(edgeintimes, newedgeintimes), pbe1$tinf.prop, pbe1$p)

  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes, newedgeintimes))
  edgeend <- c(coalnodes, edgesin, newedgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes, newedgeintimes)[nodeorder]

  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(c(edgeintimes, newedgeintimes))))[nodeorder],
                               transnode)

  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- pbe1$hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes

  ### phylotree before index case
  rewire_pullnodes_wh_strict(0)

}

rewire_pathE_wh_strict <- function() {
  ### First, dismantle minitree
  # edges entering hostID, with endtimes
  edgesin <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c")
  edgeintimes <- pbe1$v$nodetimes[edgesin]

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # coalescent nodes of hostID
  coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")

  # more coalescent nodes:
  # parent of transmission edge leaving hostID
  coalnodes <- c(take_cnode(transnode), coalnodes)

  # dismantle topology, move transmission node
  pbe1$v$nodehosts[coalnodes] <- -1L
  pbe1$v$nodehosts[transnode] <- pbe1$infector.proposed.ID
  pbe1$v$nodeparents[c(edgesin, coalnodes, transnode)] <- -1L
  pbe1$v$nodetimes[transnode] <- pbe1$tinf.prop

  ### Second, change transmission tree
  pbe1$v$inftimes[pbe1$hostID] <- pbe1$tinf.prop
  pbe1$v$infectors[pbe1$hostID] <- pbe1$infector.proposed.ID

  ### Third, rebuild minitree in hostID
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(c(edgeintimes), pbe1$tinf.prop, pbe1$p)
  coalnodes <- tail(coalnodes, -1)

  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]

  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)

  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- pbe1$hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes


  ### phylotree in infector
  rewire_pullnodes_wh_strict(pbe1$v$infectors[pbe1$hostID])
}

rewire_pathCF2_wh_strict <- function() {
  ### Identify new infector and old infector
  newinfector <- which(pbe1$v$infectors == pbe1$hostID)
  newinfector <- newinfector[which(pbe1$v$inftimes[newinfector] == min(pbe1$v$inftimes[newinfector]))]
  oldinfector <- pbe1$v$infectors[pbe1$hostID]

  ### First, remove sampling nodes and collect coalescent nodes
  # remove sample edges from new infector
  coalnodes <- c()
  sampleedges_nix <- which(pbe1$v$nodehosts == newinfector & pbe1$v$nodetypes == "x")
  for(sampleedge in sampleedges_nix) {
    coalnodes <- c(take_cnode(sampleedge), coalnodes)
  }
  coalnodes <- c(take_cnode(newinfector), coalnodes)

  # remove sample edges from hostID
  sampleedges_x <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "x")
  for(sampleedge in sampleedges_x) {
    coalnodes <- c(take_cnode(sampleedge), coalnodes)
  }
  coalnodes <- c(take_cnode(pbe1$hostID), coalnodes)

  ### Second, switch minitrees between hostID and new infector
  # transmission nodes of hostID and new infector
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
  transnode_ni <- 2*pbe1$d$nsamples - 1 + newinfector

  # switch remaining nodes in hostID and new infector
  restnodes <- which(pbe1$v$nodehosts == pbe1$hostID & (pbe1$v$nodetypes == "t" | pbe1$v$nodetypes == "c"))
  restnodes_ni <- which(pbe1$v$nodehosts == newinfector & (pbe1$v$nodetypes == "t" | pbe1$v$nodetypes == "c"))
  pbe1$v$nodehosts[restnodes] <- newinfector
  pbe1$v$nodehosts[restnodes_ni] <- pbe1$hostID

  # switch transmission nodes
  childnode_ni <- which(pbe1$v$nodeparents == transnode_ni)
  childnode <- which(pbe1$v$nodeparents == transnode)
  parentnode_ni <- pbe1$v$nodeparents[transnode_ni]
  parentnode <- pbe1$v$nodeparents[transnode]
  pbe1$v$nodeparents[childnode_ni] <- transnode
  pbe1$v$nodeparents[transnode_ni] <- parentnode
  if(parentnode_ni == transnode) {
    pbe1$v$nodeparents[transnode] <- transnode_ni
  } else {
    pbe1$v$nodeparents[childnode] <- transnode_ni
    pbe1$v$nodeparents[transnode] <- parentnode_ni
  }
  pbe1$v$nodehosts[c(transnode, transnode_ni)] <- pbe1$v$nodehosts[c(transnode_ni, transnode)]
  pbe1$v$nodetimes[c(transnode, transnode_ni)] <- pbe1$v$nodetimes[c(transnode_ni, transnode)]

  # place back sampling nodes
  pbe1$v$nodehosts[c(newinfector, sampleedges_nix)] <- newinfector
  pbe1$v$nodehosts[c(pbe1$hostID, sampleedges_x)] <- pbe1$hostID

  # Third, change transmission tree
  infectees_ni <- which(pbe1$v$infectors == newinfector)
  infectees <- which(pbe1$v$infectors == pbe1$hostID)
  pbe1$v$inftimes[c(pbe1$hostID, newinfector)] <- pbe1$v$inftimes[c(newinfector, pbe1$hostID)]
  pbe1$v$infectors[infectees] <- newinfector
  pbe1$v$infectors[infectees_ni] <- pbe1$hostID
  pbe1$v$infectors[c(pbe1$hostID, newinfector)] <- c(newinfector, oldinfector)

  # Fourth, add sample edges in hostID and new infector
  rewire_pullnodes_wh_strict(pbe1$hostID)
  rewire_pullnodes_wh_strict(newinfector)
  rewire_pullnodes_wh_strict(oldinfector)
}

rewire_pathK_wh_strict <- function() {
  ### First, dismantle minitree
  # edges entering hostID, with endtimes
  edgesin <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c")
  edgeintimes <- pbe1$v$nodetimes[edgesin]

  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID

  # all coalescent nodes in new infector and hostID
  coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")

  # dismantle topology, move transmission node
  pbe1$v$nodehosts[coalnodes] <- -1
  pbe1$v$nodeparents[c(edgesin, coalnodes)] <- -1

  ### Second, rebuild minitree
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(edgeintimes, pbe1$v$inftimes[pbe1$hostID], pbe1$p)

  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]

  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)

  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- pbe1$hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes
}

rewire_pullnodes_wh_strict <- function(currentID) {
  loosenodes <- which(pbe1$v$nodehosts == currentID & pbe1$v$nodeparents == -1)
  if(length(loosenodes) > 0) {
    free_cnodes <- which(pbe1$v$nodetypes == "c" & pbe1$v$nodeparents == -1)
    if(currentID == 0) {
      pbe1$v$nodeparents[loosenodes] <- 0L
     } else {
      # old edges entering currentID, with endtimes
      edgesendold <- setdiff(which(pbe1$v$nodehosts == currentID & pbe1$v$nodetypes != "c"), loosenodes)
      if(length(edgesendold) == 0) {
        parentnode <- 2 * pbe1$d$nsamples - 1 + currentID
        edgesendold <- currentID
        pbe1$v$nodeparents[edgesendold] <- parentnode
        pbe1$v$nodehosts[parentnode] <- pbe1$v$infectors[currentID]
        pbe1$v$nodetimes[parentnode] <- pbe1$v$inftimes[currentID]
        loosenodes <- setdiff(loosenodes, currentID)
      }
      edgesendoldtimes <- pbe1$v$nodetimes[edgesendold]

      if(length(loosenodes) > 0) {
        # coalescentnodes in currentID, with endtimes
        coalescentnodesold <- which(pbe1$v$nodehosts == currentID & pbe1$v$nodetypes == "c")
        coalescenttimesold <- pbe1$v$nodetimes[coalescentnodesold]

        # endtimes of new edges, and new coalescenttimes
        loosenodetimes <- c()
        coalescenttimesnew <- c()
        for(le in 1:length(loosenodes)) {
          coalescenttimesnew <- c(coalescenttimesnew,
                                    sample_singlecoaltime(c(edgesendoldtimes, loosenodetimes),
                                                          c(coalescenttimesold, coalescenttimesnew),
                                                          pbe1$v$nodetimes[loosenodes[le]], pbe1$v$inftimes[currentID], pbe1$p))
          loosenodetimes <- c(loosenodetimes, pbe1$v$nodetimes[loosenodes[le]])
        }

        # old within-host minitree
        childnodes <- c(edgesendold, coalescentnodesold)
        parentnodes <- pbe1$v$nodeparents[childnodes]
        childnodestimes <- c(edgesendoldtimes, coalescenttimesold)

        # place new edges into minitree
        loosenodestoinfector <- c()
        free_cnodestoinfector <- c()
        for(le in 1:length(loosenodes)) {
            newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
            childnodes <- c(childnodes, loosenodes[le], free_cnodes[le])
            parentnodes <- c(parentnodes, free_cnodes[le], parentnodes[childnodes == newchildnode])
            childnodestimes <- c(childnodestimes, loosenodetimes[le], coalescenttimesnew[le])
            parentnodes[childnodes == newchildnode] <- free_cnodes[le]
        }

        # change phybreak object
        pbe1$v$nodetimes[childnodes] <- childnodestimes
        pbe1$v$nodehosts[childnodes] <- currentID
        pbe1$v$nodeparents[childnodes] <- parentnodes
      }
    }
  }
}
