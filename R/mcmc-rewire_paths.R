rewire_pathA <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathA_wh_loose()
  } else {
    rewire_pathA_wh_strict()
  }
}

rewire_pathB <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathB_wh_loose()
  } else {
    rewire_pathB_wh_strict()
  }
}

rewire_pathCF1 <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathCF1_wh_loose()
  } else {
    rewire_pathCF1_wh_strict()
  }
}

rewire_pathD <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathD_wh_loose()
  } else {
    rewire_pathD_wh_strict()
  }
}

rewire_pathE <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathE_wh_loose()
  } else {
    rewire_pathE_wh_strict()
  }
}

rewire_pathCF2 <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathCF2_wh_loose()
  } else {
    rewire_pathCF2_wh_strict()
  }
}

rewire_pathK <- function(loose_bottleneck) {
  if(loose_bottleneck) {
    rewire_pathK_wh_loose()
  } else {
    rewire_pathK_wh_strict()
  }
}


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

take_cnode <- function(childnode) {
  parentnode <- pbe1$v$nodeparents[childnode]
  while(pbe1$v$nodetypes[parentnode] %in% c("t", "b")) {
    pbe1$v$nodehosts[childnode] <- -1
    pbe1$v$nodeparents[childnode] <- -1
    
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
