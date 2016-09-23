### functions to make objects of class phylo ###


### make a class 'phylo' object from class 'phybreak' input called by: get.phylo sim.phybreak required packages: phangorn
.makephylo.phybreak <- function(nodetimes, nodeparents, nodenames) {
    ### topology
    Nhosts <- (1 + length(nodetimes))/3
    indext <- (1:length(nodetimes))[nodeparents == 0]
    indexc <- (1:length(nodetimes))[nodeparents == indext]
    edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
    edgeend <- (1:length(nodetimes))[nodeparents != 0 & nodeparents != indext]
    edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]
    
    toremove <- (edgeend >= 2 * Nhosts)
    toreplace <- (edgestart >= 2 * Nhosts)
    edgelengths[toreplace] <- edgelengths[toreplace] + edgelengths[toremove][match(edgestart[toreplace], edgeend[toremove])]
    edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace], edgeend[toremove])]
    edgestart <- edgestart[!toremove]
    edgeend <- edgeend[!toremove]
    edgelengths <- edgelengths[!toremove]
    
    ### give first coalescent node number Nhosts+1
    if (indexc != Nhosts + 1) {
        edgestart[edgestart == indexc] <- 0
        edgeend[edgeend == indexc] <- 0
        edgestart[edgestart == Nhosts + 1] <- indexc
        edgeend[edgeend == Nhosts + 1] <- indexc
        edgestart[edgestart == 0] <- Nhosts + 1
        edgeend[edgeend == 0] <- Nhosts + 1
    }
    
    edges <- matrix(c(edgestart, edgeend), ncol = 2)
    
    toname <- nodenames
    
    
    res <- list(edge = edges, edge.length = edgelengths, Nnode = Nhosts - 1, tip.label = toname)
    class(res) <- "phylo"
    res <- ape::reorder.phylo(res)
    res <- ape::ladderize(res)
    return(res)
    
}


### make a class 'phylo' object from class 'phybreak' input, including Simmap items called by: get.phylo required packages:
### phytools phangorn
.makephylosimmap.phybreak <- function(nodetimes, nodeparents, nodehosts, nodenames) {
    ### topology for phylo
    Nhosts <- (1 + length(nodetimes))/3
    indext <- (1:length(nodetimes))[nodeparents == 0]
    indexc <- (1:length(nodetimes))[nodeparents == indext]
    edgestart <- nodeparents[nodeparents != 0 & nodeparents != indext]
    edgeend <- (1:length(nodetimes))[nodeparents != 0 & nodeparents != indext]
    edgelengths <- nodetimes[edgeend] - nodetimes[edgestart]
    
    toremove <- (edgeend >= 2 * Nhosts)
    toreplace <- (edgestart >= 2 * Nhosts)
    edgelengths[toreplace] <- edgelengths[toreplace] + edgelengths[toremove][match(edgestart[toreplace], edgeend[toremove])]
    edgestart[toreplace] <- edgestart[toremove][match(edgestart[toreplace], edgeend[toremove])]
    edgestart <- edgestart[!toremove]
    edgeend <- edgeend[!toremove]
    edgelengths <- edgelengths[!toremove]
    
    edges <- matrix(c(edgestart, edgeend), ncol = 2)
    
    ### items for simmap (phytools)
    whichcolor <- function(iors, path) {
        if (path[1] == 0) {
            return(length(path)%%3)
        } else {
            return(whichcolor(iors, c(iors[path[1]], path)))
        }
    }
    
    tipstates <- c("red", "blue", "black")[1 + sapply(1:Nhosts, whichcolor, iors = tail(nodehosts, Nhosts))]
    nodestates <- matrix(tipstates[nodehosts[edges]], ncol = 2)
    edgemaps <- list()
    for (i in 1:(2 * (Nhosts - 1))) {
        if (nodestates[i, 1] == nodestates[i, 2]) {
            edgemaps[[i]] <- edgelengths[i]
            names(edgemaps[[i]]) <- nodestates[i, 1]
        } else {
            nodes <- edges[i, ]
            ntimes <- nodetimes[c(nodes[1], nodeparents[nodes[2]], nodes[2])]
            edgemaps[[i]] <- ntimes[2:3] - ntimes[1:2]
            names(edgemaps[[i]]) <- nodestates[i, ]
        }
    }
    
    mappededge <- matrix(nrow = 2 * (Nhosts - 1), ncol = 3)
    for (i in 1:(2 * (Nhosts - 1))) {
        mappededge[i, ] <- edgemaps[[i]][c("blue", "black", "red")]
    }
    mappededge[is.na(mappededge)] <- 0
    colnames(mappededge) <- c("blue", "black", "red")
    rownames(mappededge) <- paste0(edges[, 1], ",", edges[, 2])
    
    
    ### give first coalescent node number Nhosts+1
    if (indexc != Nhosts + 1) {
        edgestart[edgestart == indexc] <- 0
        edgeend[edgeend == indexc] <- 0
        edgestart[edgestart == Nhosts + 1] <- indexc
        edgeend[edgeend == Nhosts + 1] <- indexc
        edgestart[edgestart == 0] <- Nhosts + 1
        edgeend[edgeend == 0] <- Nhosts + 1
    }
    edges <- matrix(c(edgestart, edgeend), ncol = 2)
    
    ### make preliminary result for ladderizing
    respreladder <- list(edge = edges, edge.length = edgelengths, Nnode = Nhosts - 1, tip.label = nodenames)
    class(respreladder) <- c("phylo")
    respreladder <- ape::reorder.phylo(respreladder)
    
    ### ladderize to determine reordering of edges
    edgeladder <- ape::ladderize(respreladder)$edge[, 2]
    newedgeorder <- match(edgeladder, edgeend)
    newtiporder <- match(respreladder$tip.label, nodenames)
    
    reorderedmaps <- list()
    for (i in 1:(2 * (Nhosts - 1))) {
        reorderedmaps[[i]] <- edgemaps[[newedgeorder[i]]]
    }
    
    
    res <- list(edge = edges[newedgeorder, ], edge.length = edgelengths[newedgeorder], Nnode = Nhosts - 1, tip.label = nodenames[newtiporder], 
        node.state = nodestates[newedgeorder, ], states = tipstates[newtiporder], maps = reorderedmaps, mapped.edge = mappededge[newedgeorder, 
            ])
    class(res) <- c("simmap", "phylo")
    attr(res, "order") <- "cladewise"
    
    return(res)
    
}
