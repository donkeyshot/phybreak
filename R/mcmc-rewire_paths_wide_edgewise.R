rewire_change_infector_wide <- function(ID, newinfector) {
  btnodes <- pbe1$v$nodeparents[pbe1$v$nodehosts == ID]
  btnodes <- btnodes[btnodes > 0]
  btnodes <- btnodes[pbe1$v$nodetypes[btnodes] != "c"]
  
  if(pbe1$v$infectors[ID] == 0) {
    coalnodes <- which(pbe1$v$nodehosts == 0 & pbe1$v$nodetypes == "c")
    pbe1$v$nodeparents[c(btnodes, coalnodes)] <- -1L
    pbe1$v$nodehosts[coalnodes] <- -1
  } else {
    coalnodes <- sapply(btnodes, take_cnode)
  }

  pbe1$v$nodehosts[btnodes] <- newinfector
  pbe1$v$infectors[ID] <- newinfector
  return(coalnodes)
}


rewire_within_wide_edgewise <- function(ID, tinf) {
  coalnodes <- which(pbe1$v$nodehosts == ID & pbe1$v$nodetypes == "c")
  coaltimes_old <- pbe1$v$nodetimes[coalnodes]
  
  coaltimes_new <- sample_coaltimes(pbe1$v$nodetimes[pbe1$v$nodehosts == ID & pbe1$v$nodetypes != "c"],
                                    tinf, pbe1$p)
  coaltimes_new <- coaltimes_new[coaltimes_new > tinf]
  

  if(length(coaltimes_new) == length(coaltimes_old)) {
    pbe1$logLiktoporatio <- pbe1$logLiktoporatio - lik_topology_host(pbe1, ID)
    
    coaltimes_new <- coaltimes_new[rank(coaltimes_old)]
    
    pbe1$v$nodetimes[coalnodes] <- coaltimes_new
    btnodes <- pbe1$v$nodeparents[pbe1$v$nodehosts == ID]
    btnodes <- btnodes[pbe1$v$nodetypes[btnodes] != "c"]
    pbe1$v$nodetimes[btnodes] <- tinf
    pbe1$v$inftimes[ID] <- tinf
    
    if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == ID] - 
           pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == ID]] < 0)) {
      pbe1$logLiktoporatio <- -Inf
      return()
    } 
    
    pbe1$logLiktoporatio <- pbe1$logLiktoporatio + lik_topology_host(pbe1, ID)
  } else if(length(coaltimes_new) < length(coaltimes_old)) {
    cnodediff <- length(coaltimes_old) - length(coaltimes_new)
    toremove <- coalnodes[rank(coaltimes_old) <= cnodediff]
    newrootnodes <- setdiff(which(pbe1$v$nodeparents %in% toremove), toremove)
    rootparents <- c(setdiff(pbe1$v$nodeparents[toremove], toremove), 
                     replicate(cnodediff, rewire_selectbnode(pbe1$v$infectors[ID])))
    if(any(newrootnodes %in% .ptr(pbe1$v$nodeparents, ID))) {
      rootparents <- link_s_to_t(rootparents, newrootnodes, 
                                 newrootnodes[newrootnodes %in% .ptr(pbe1$v$nodeparents, ID)], 
                                 ID + 2*pbe1$d$nsamples - 1)
    }
    pbe1$v$nodeparents[newrootnodes] <- rootparents
    pbe1$v$nodeparents[toremove] <- -1L
    pbe1$v$nodehosts[toremove] <- -1L
    
    coalnodes <- which(pbe1$v$nodehosts == ID & pbe1$v$nodetypes == "c")
    coaltimes_old <- pbe1$v$nodetimes[coalnodes]
    
    pbe1$logLiktoporatio <- pbe1$logLiktoporatio - lik_topology_host(pbe1, ID)
    
    coaltimes_new <- coaltimes_new[rank(coaltimes_old)]
    
    pbe1$v$nodetimes[coalnodes] <- coaltimes_new
    btnodes <- pbe1$v$nodeparents[pbe1$v$nodehosts == ID]
    btnodes <- btnodes[pbe1$v$nodetypes[btnodes] != "c"]
    pbe1$v$nodetimes[btnodes] <- tinf
    pbe1$v$inftimes[ID] <- tinf
    
    if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == ID] - 
           pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == ID]] < 0)) {
      pbe1$logLiktoporatio <- -Inf
      return()
    } 
    
    pbe1$logLiktoporatio <- pbe1$logLiktoporatio + lik_topology_host(pbe1, ID)
  } else {
    pbe1$logLiktoporatio <- pbe1$logLiktoporatio - lik_topology_host(pbe1, ID)
    
    cnodediff <- length(coaltimes_new) - length(coaltimes_old)
    coaltimes_toadd <- coaltimes_new[1:cnodediff]
    coaltimes_new <- tail(coaltimes_new, -cnodediff)[rank(coaltimes_old)]
    
    pbe1$v$nodetimes[coalnodes] <- coaltimes_new
    btnodes <- pbe1$v$nodeparents[pbe1$v$nodehosts == ID]
    btnodes <- sort(btnodes[pbe1$v$nodetypes[btnodes] != "c"])
    pbe1$v$nodetimes[btnodes] <- tinf
    pbe1$v$inftimes[ID] <- tinf
    
    if(any(pbe1$v$nodetimes[pbe1$v$nodehosts == ID] - 
           pbe1$v$nodetimes[pbe1$v$nodeparents[pbe1$v$nodehosts == ID]] < 0)) {
      pbe1$logLiktoporatio <- -Inf
      return()
    } 
    
    pbe1$logLiktoporatio <- pbe1$logLiktoporatio + lik_topology_host(pbe1, ID)
    
    toadd <- tail(which(pbe1$v$nodetypes == "c"& pbe1$v$nodehosts == -1), cnodediff)
    pbe1$v$nodetimes[toadd] <- coaltimes_toadd
    oldrootnodes <- which(pbe1$v$nodehosts == ID)
    oldrootnodes <- oldrootnodes[pbe1$v$nodetypes[pbe1$v$nodeparents[oldrootnodes]] != "c"]
    rootparents <- sort(pbe1$v$nodeparents[oldrootnodes])
    for(i in tail(rootparents, cnodediff)) rewire_removebnode(i)
    rootparents <- head(rootparents, -cnodediff)
    nodes_newtopo <- c(toadd, oldrootnodes)
    nodetimes_newtopo <- pbe1$v$nodetimes[nodes_newtopo]
    nodetypes_newtopo <- c(rep("c", cnodediff), rep("x", length(oldrootnodes)))
    nodes_newtopo <- nodes_newtopo[order(nodetimes_newtopo)]
    nodetypes_newtopo <- nodetypes_newtopo[order(nodetimes_newtopo)]
    nodetimes_newtopo <- sort(nodetimes_newtopo)
    nodeparents_newtopo <- sample_topology(nodes_newtopo, nodetimes_newtopo, nodetypes_newtopo, rootparents)
    nodeparents_newtopo <- link_s_to_t(nodeparents_newtopo, nodes_newtopo, 
                                       nodes_newtopo[nodes_newtopo %in% .ptr(pbe1$v$nodeparents, ID)], 
                                       ID + 2*pbe1$d$nsamples - 1)
    pbe1$v$nodeparents[nodes_newtopo] <- nodeparents_newtopo
    pbe1$v$nodehosts[nodes_newtopo] <- ID
  }
}

rewire_selectbnode <- function(infector) {
  if(any(pbe1$v$nodetypes == "0")) {
    bnode <- which(pbe1$v$nodetypes == "0")[1]
    pbe1$v$nodetypes[bnode] <- "b"
    pbe1$v$nodehosts[bnode] <- infector
  } else {
    bnode <- length(pbe1$v$nodetypes) + 1L
    pbe1$v$nodetypes <- c(pbe1$v$nodetypes, "b")
    pbe1$v$nodehosts <- c(pbe1$v$nodehosts, infector)
    pbe1$v$nodeparents <- c(pbe1$v$nodeparents, -1L)
    pbe1$v$nodetimes <- c(pbe1$v$nodetimes, 0)
  }
  return(bnode)
}

rewire_removebnode <- function(bnode) {
  pbe1$v$nodetypes[bnode] <- "0"
  pbe1$v$nodehosts[bnode] <- -1L
  pbe1$v$nodeparents[bnode] <- -1L
}

rewire_pathA_wide_edgewise <- function() {
  rewire_change_infector_wide(pbe1$hostID, 0L)
  rewire_within_wide_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_pullnodes_wide(0L)
  }
}


rewire_pathB_wide_edgewise <- function() {
  ### Identify new index
  newindex <- which(pbe1$v$inftimes == sort(pbe1$v$inftimes)[2])

  rewire_change_infector_wide(newindex, 0L)
  rewire_change_infector_wide(pbe1$hostID, pbe1$infector.proposed.ID)
  rewire_within_wide_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_pullnodes_wide(0)
    rewire_pullnodes_wide(pbe1$infector.proposed.ID)
  }
}

rewire_pathCF_wide_edgewise <- function() {
  ### Identify new infector and old infector
  newinfector <- which(pbe1$v$infectors == pbe1$hostID)
  newinfector <- newinfector[which(pbe1$v$inftimes[newinfector] == min(pbe1$v$inftimes[newinfector]))]
  oldinfector <- pbe1$v$infectors[pbe1$hostID]
  
  firstinfectiontime <- pbe1$v$inftimes[pbe1$hostID]
  secondinfectiontime <- pbe1$v$inftimes[newinfector]
  
  rewire_change_infector_wide(newinfector, oldinfector)
  rewire_within_wide_edgewise(newinfector, firstinfectiontime)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_change_infector_wide(pbe1$hostID, newinfector)
    rewire_within_wide_edgewise(pbe1$hostID, secondinfectiontime)
  }
  if(pbe1$logLiktoporatio > -Inf) rewire_pullnodes_wide(newinfector)
}

rewire_pathD_wide_edgewise <- function() {
  currentindex <- which(pbe1$v$infectors == 0)
  
  rewire_change_infector_wide(pbe1$hostID, 0L)
  rewire_within_wide_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) {
    rewire_change_infector_wide(currentindex, pbe1$hostID)
    rewire_pullnodes_wide(pbe1$hostID)
  }
}

rewire_pathE_wide_edgewise <- function() {
  rewire_change_infector_wide(pbe1$hostID, pbe1$infector.proposed.ID)
  rewire_within_wide_edgewise(pbe1$hostID, pbe1$tinf.prop)
  if(pbe1$logLiktoporatio > -Inf) rewire_pullnodes_wide(pbe1$infector.proposed.ID)
}
  

rewire_pathK_wide_classic <- function() {
  ### First, dismantle minitree
  # edges entering hostID, with endtimes
  edgesin <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c")
  edgeintimes <- pbe1$v$nodetimes[edgesin]
  
  # transmission node of hostID
  transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
  
  # bottleneck nodes
  all_edges <- which(pbe1$v$nodehosts == pbe1$hostID)
  bottlenodes <- pbe1$v$nodeparents[all_edges]
  bottlenodes <- bottlenodes[pbe1$v$nodetypes[bottlenodes] == "b"]
  
  # all coalescent nodes in new infector and hostID
  coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
  
  # more coalescent nodes:
  # parents of bottleneck edges leaving hostID
  for(bnode in bottlenodes) {
    coalnodes <- c(take_cnode(bnode), coalnodes)
  }
  
  # dismantle topology, move transmission node
  pbe1$v$nodetypes[bottlenodes] <- 0
  pbe1$v$nodehosts[c(coalnodes, bottlenodes)] <- -1
  pbe1$v$nodeparents[c(edgesin, coalnodes, bottlenodes)] <- -1
  
  ### Second, rebuild minitree
  # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
  newcoaltimes <- sample_coaltimes(edgeintimes, pbe1$v$inftimes[pbe1$hostID], pbe1$p)
  nr_of_bnodes <- sum(newcoaltimes < pbe1$v$inftimes[pbe1$hostID])
  if(nr_of_bnodes > 0) {
    newcoaltimes <- tail(newcoaltimes, - nr_of_bnodes)
  }
  coalnodes <- tail(coalnodes, length(newcoaltimes))
  
  # order all edges (transmission, sample, coalescent) by endtime in hostID
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]
  
  # take bottleneck nodes as needed
  if(nr_of_bnodes > 0) {
    availablebnodes <- which(pbe1$v$nodetypes == "0")
    if(length(availablebnodes) >= nr_of_bnodes) {
      bnodes <- head(availablebnodes, nr_of_bnodes)
    } else {
      bnode_shortage <- nr_of_bnodes - length(availablebnodes)
      bnodes <- c(availablebnodes, length(pbe1$v$nodetypes) + 1:bnode_shortage)
    }
    pbe1$v$nodetypes[bnodes] <- "b"
    pbe1$v$nodehosts[bnodes] <- pbe1$v$infectors[pbe1$hostID]
    pbe1$v$nodeparents[bnodes] <- -1L
    pbe1$v$nodetimes[bnodes] <- pbe1$v$inftimes[pbe1$hostID]
  } else {
    bnodes <- c()
  }
  
  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend, 
                               edgeendtimes, 
                               c(rep("c", length(newcoaltimes)), 
                                 rep("x", length(edgeintimes)))[nodeorder],
                               c(transnode, bnodes))
  if(nr_of_bnodes > 0) edgestart <- link_s_to_t(edgestart, edgeend, 
                                                snode = pbe1$hostID, 
                                                tnode = transnode)
  
  # change minitree in hostID
  pbe1$v$nodehosts[edgeend] <- pbe1$hostID
  pbe1$v$nodeparents[edgeend] <- edgestart
  pbe1$v$nodetimes[edgeend] <- edgeendtimes
  
  ### phylotree in infector
  rewire_pullnodes_wide(pbe1$v$infectors[pbe1$hostID])
}




