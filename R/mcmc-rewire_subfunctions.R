### attach ("pull") new tips to the existing minitree in currentID (complete bottleneck)
rewire_pullnodes_complete <- function(currentID) {
  
  # tips to be connected in currentID
  loosenodes <- which(pbe1$v$nodehosts == currentID & pbe1$v$nodeparents == -1)
  
  # only if there is a loose node
  if(length(loosenodes) > 0) {
    
    # available coalescent nodes
    free_cnodes <- which(pbe1$v$nodetypes == "c" & pbe1$v$nodeparents == -1)
    
    # if(currentID == 0 && length(free_cnodes) == 0) {
    #   if(sum(pbe1$v$infectors==0) == 1){
    #     pbe1$v$nodeparents[loosenodes] <- 0L
    #   } else {
    #     coalnodes <- which(pbe1$v$nodetypes == "c" & pbe1$v$nodehosts == 0)
    #     free_cnodes <- coalnodes[sapply(coalnodes, function(x) sum(pbe1$v$nodeparents == x) == 1)]
    #     pbe1$v$nodeparents[loosenodes] <- free_cnodes
    #   }
    #   return()
    # }
    
    # edges already entering currentID
    edgesendold <- setdiff(which(pbe1$v$nodehosts == currentID & pbe1$v$nodetypes != "c"), loosenodes)
    
    # if no edges already entering, make one by connecting one tip to host's root node
    if(length(edgesendold) == 0) {
      # in history host: connect first tip to non-existing node 0
      if(currentID == 0) {
        parentnode <- 0
        edgesendold <- loosenodes[1]
        pbe1$v$nodeparents[edgesendold] <- parentnode
        loosenodes <- setdiff(loosenodes, edgesendold)
      } else 
        # in other hosts, connect sampling tip to infection node
      {
        parentnode <- 2 * pbe1$d$nsamples - 1 + currentID
        edgesendold <- currentID
        pbe1$v$nodeparents[edgesendold] <- parentnode
        pbe1$v$nodehosts[parentnode] <- pbe1$v$infectors[currentID]
        pbe1$v$nodetimes[parentnode] <- pbe1$v$inftimes[currentID]
        loosenodes <- setdiff(loosenodes, edgesendold)
      }
    }
    
    # old edges entering currentID, with endtimes
    edgesendoldtimes <- pbe1$v$nodetimes[edgesendold]
    
    # only if any loose node still present
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
                                                      pbe1$v$nodetimes[loosenodes[le]], pbe1$v$inftimes[currentID], pbe1$p, currentID == 0))
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

### attach ("pull") new tips to the existing minitree in currentID (wide bottleneck)
rewire_pullnodes_wide <- function(currentID) {
  loosenodes <- which(pbe1$v$nodehosts == currentID & pbe1$v$nodeparents == -1)
  if(length(loosenodes) > 0) {
    free_cnodes <- which(pbe1$v$nodetypes == "c" & pbe1$v$nodeparents == -1)
    if(currentID == 0) {
      # identify index case
      indexID <- which(pbe1$v$infectors == 0)
      
      # old edges entering currentID, with endtimes
      edgesendold <- setdiff(which(pbe1$v$nodehosts == currentID & (pbe1$v$nodetypes == "t" | pbe1$v$nodetypes == "b")), 
                             loosenodes)
      if(length(edgesendold) == 0) {
        edgesendold <- 2 * pbe1$d$nsamples - 1 + indexID
        pbe1$v$nodeparents[edgesendold] <- 0L
        loosenodes <- setdiff(loosenodes, edgesendold)
      }
      edgesendoldtimes <- pbe1$v$nodetimes[edgesendold]
      
      if(length(loosenodes) > 0) {
        # coalescentnodes in currentID, with endtimes
        coalescentnodesold <- setdiff(which(pbe1$v$nodehosts == currentID & pbe1$v$nodetypes == "c"), free_cnodes)
        coalescenttimesold <- pbe1$v$nodetimes[coalescentnodesold]
        
        # endtimes of new edges, and new coalescenttimes
        loosenodetimes <- c()
        coalescenttimesnew <- c()
        for(le in 1:length(loosenodes)) {
          coalescenttimesnew <- c(coalescenttimesnew, 
                                  sample_singlecoaltime(c(edgesendoldtimes, loosenodetimes),
                                                        c(coalescenttimesold, coalescenttimesnew),
                                                        pbe1$v$nodetimes[loosenodes[le]], pbe1$v$inftimes[indexID] - pbe1$p$sample.mean, pbe1$p))
          loosenodetimes <- c(loosenodetimes, pbe1$v$nodetimes[loosenodes[le]])
        }
        
        # old within-host minitree
        childnodes <- c(edgesendold, coalescentnodesold)
        parentnodes <- pbe1$v$nodeparents[childnodes]
        childnodestimes <- c(edgesendoldtimes, coalescenttimesold)
        
        # place new edges into minitree
        for(le in 1:length(loosenodes)) {
          newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
          childnodes <- c(childnodes, loosenodes[le], free_cnodes[le])
          parentnodes <- c(parentnodes, free_cnodes[le], parentnodes[childnodes == newchildnode])
          childnodestimes <- c(childnodestimes, loosenodetimes[le], coalescenttimesnew[le])
          parentnodes[childnodes == newchildnode] <- free_cnodes[le]
        }
        
        # change phybreak object
        pbe1$v$nodetimes[childnodes] <- childnodestimes
        pbe1$v$nodehosts[childnodes] <- 0L
        pbe1$v$nodeparents[childnodes] <- parentnodes
      }
      
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
        
        # take bottleneck nodes as needed
        if((nr_of_bnodes <- sum(coalescenttimesnew < pbe1$v$inftimes[currentID])) > 0) {
          availablebnodes <- which(pbe1$v$nodetypes == "0")
          availablebnodes <- availablebnodes[availablebnodes >= 2 * pbe1$d$nsamples + pbe1$p$obs]
          if(length(availablebnodes) >= nr_of_bnodes) {
            bnodes <- head(availablebnodes, nr_of_bnodes)
          } else {
            bnode_shortage <- nr_of_bnodes - length(availablebnodes)
            bnodes <- c(availablebnodes, length(pbe1$v$nodeparents) + 1:bnode_shortage)
            pbe1$v$nodetypes <- c(pbe1$v$nodetypes, rep("b", bnode_shortage)) # other nodevectors are extended automatically
          }
          pbe1$v$nodehosts[bnodes] <- pbe1$v$infectors[currentID]
          pbe1$v$nodetimes[bnodes] <- pbe1$v$inftimes[currentID]
          pbe1$v$nodetypes[bnodes] <- "b"
          pbe1$v$nodeparents[bnodes] <- -1
        }
        
        # place new edges into minitree
        loosenodestoinfector <- c()
        free_cnodestoinfector <- c()
        for(le in 1:length(loosenodes)) {
          if(coalescenttimesnew[le] >= pbe1$v$inftimes[currentID]) {
            newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
            childnodes <- c(childnodes, loosenodes[le], free_cnodes[le])
            parentnodes <- c(parentnodes, free_cnodes[le], parentnodes[childnodes == newchildnode])
            childnodestimes <- c(childnodestimes, loosenodetimes[le], coalescenttimesnew[le])
            parentnodes[childnodes == newchildnode] <- free_cnodes[le]
          } else {
            childnodes <- c(childnodes, loosenodes[le])
            parentnodes <- c(parentnodes, bnodes[1])
            childnodestimes <- c(childnodestimes, loosenodetimes[le])
            bnodes <- bnodes[-1]
          }
        }
        parentnodes <- link_s_to_t(parentnodes, childnodes, currentID, 2 * pbe1$d$nsamples - 1 + currentID)
        
        # change phybreak object
        pbe1$v$nodetimes[childnodes] <- childnodestimes
        pbe1$v$nodehosts[childnodes] <- currentID
        pbe1$v$nodeparents[childnodes] <- parentnodes
        
        rewire_pullnodes_wide(pbe1$v$infectors[currentID])
        
      }
    }
  }
}

### remove tips from existing minitree, and return the coalescent node
take_cnode <- function(childnode) {
  parentnode <- pbe1$v$nodeparents[childnode]
  while(pbe1$v$nodetypes[parentnode] %in% c("t", "b")) {
    pbe1$v$nodehosts[childnode] <- -1L
    pbe1$v$nodeparents[childnode] <- -1L
    
    childnode <- parentnode
    parentnode <- pbe1$v$nodeparents[childnode]
    if(pbe1$v$nodetypes[childnode] == "b") pbe1$v$nodetypes[childnode] <- "0"
  }
  second_childnode <- setdiff(which(pbe1$v$nodeparents == parentnode), childnode)

  pbe1$v$nodeparents[second_childnode] <- pbe1$v$nodeparents[parentnode]
  pbe1$v$nodehosts[c(childnode, parentnode)] <- -1L
  pbe1$v$nodeparents[c(childnode, parentnode)] <- -1L
  return(parentnode)
}

### switch t and b nodes to let the s-node have the t-node as a direct ancestor
link_s_to_t <- function(parentnodes, childnodes, snode, tnode) {
  current_t_node <- which(parentnodes == tnode)
  
  find_current_tbnode <- which(snode == childnodes)
  while(parentnodes[find_current_tbnode] %in% childnodes) {
    find_current_tbnode <- which(parentnodes[find_current_tbnode] == childnodes)
  }
  
  if(current_t_node != find_current_tbnode) {
    parentnodes[current_t_node] <- parentnodes[find_current_tbnode]
    parentnodes[find_current_tbnode] <- tnode
  }
  
  return(parentnodes)
}

rearrange_tree_matrix <- function(){
  v <- pbe1$v
  p <- pbe1$p

  v$tree <- sapply(1:p$obs, function(x) tail(.ptr(v$infectors, x), 1))
  
  copy2pbe1("v", environment())
}
