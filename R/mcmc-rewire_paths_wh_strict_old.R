# rewire_pathA_wh_strict <- function() {
#   ### Make input locally available
#   d <- pbe1$d
#   p <- pbe1$p
#   v <- pbe1$v
#   hostID <- pbe1$hostID
#   tinf.prop <- pbe1$tinf.prop
# 
#   ### Bookkeeping: change transmission tree
#   v$nodetimes[2*d$nsamples - 1 + hostID] <- tinf.prop
#   v$inftimes[hostID] <- tinf.prop
# 
# 
#   ### Change minitree
#   # edges entering hostID, with endtimes
#   edgesin <- which(v$nodehosts == hostID & v$nodetypes != "c")
#   edgeintimes <- v$nodetimes[edgesin]
# 
#   # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
#   v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-
#     sample_coaltimes(edgeintimes, tinf.prop, p)
# 
#   # order all edges (transmission, sample, coalescent) by endtime in hostID
#   edgeendtimes <- v$nodetimes[v$nodehosts == hostID]
#   edgeendtypes <- v$nodetypes[v$nodehosts == hostID]
#   nodeorder <- order(edgeendtimes, edgeendtypes)
#   edgesend <- which(v$nodehosts == hostID)[nodeorder]
# 
#   # sample topology of minitree within hostID
#   v$nodeparents[edgesend] <-
#     sample_topology(edgesend,
#                     edgeendtimes[nodeorder],
#                     edgeendtypes[nodeorder],
#                     2*d$nsamples - 1 + hostID)
# 
#   copy2pbe1("v", environment())
# }
# 
# rewire_pathB_wh_strict <- function() {
#   ### Make input locally available
#   d <- pbe1$d
#   p <- pbe1$p
#   v <- pbe1$v
#   hostID <- pbe1$hostID
#   tinf.prop <- pbe1$tinf.prop
#   infector.proposed.ID <- pbe1$infector.proposed.ID
# 
#   ### Identify new index
#   newindex <- which(v$inftimes == sort(v$inftimes)[2])
#   transnode_ni <- 2*d$nsamples - 1 + newindex
#   transnode <- 2*d$nsamples - 1 + hostID
# 
#   ### Bookkeeping: change infection time
#   v$nodetimes[transnode] <- tinf.prop
#   v$inftimes[hostID] <- tinf.prop
# 
#   ### Move coalescent node from hostID to new infector
#   # remove from hostID
#   movingcoalnode <- v$nodeparents[transnode_ni]
#   otherchildnode <- setdiff(which(v$nodeparents == movingcoalnode), transnode_ni)
#   v$nodeparents[c(transnode_ni, transnode, otherchildnode)] <-
#     c(0L, movingcoalnode, v$nodeparents[movingcoalnode])
#   # sample time and edge (child node) in new infector
#   v$nodetimes[movingcoalnode] <-
#     sample_singlecoaltime(v$nodetimes[v$nodehosts == infector.proposed.ID & v$nodetypes != "c"],
#                           v$nodetimes[v$nodehosts == infector.proposed.ID & v$nodetypes == "c"],
#                           v$inftimes[hostID],
#                           v$inftimes[infector.proposed.ID],
#                           p)
#   otherchildnode <-
#     sample_singlechildnode(which(v$nodehosts == infector.proposed.ID),
#                            v$nodeparents[v$nodehosts == infector.proposed.ID],
#                            v$nodetimes[v$nodehosts == infector.proposed.ID],
#                            v$nodetimes[movingcoalnode])
#   # place in new infector
#   v$nodeparents[c(movingcoalnode, otherchildnode)] <-
#     c(v$nodeparents[otherchildnode], movingcoalnode)
#   v$nodehosts[movingcoalnode] <- infector.proposed.ID
# 
#   ### Bookkeeping: change transmission tree topology
#   v$nodehosts[c(transnode_ni, transnode)] <- c(0L, infector.proposed.ID)
#   v$infectors[c(newindex, hostID)] <- c(0L, infector.proposed.ID)
# 
#   ### Change minitree
#   # edges entering hostID, with endtimes
#   edgesin <- which(v$nodehosts == hostID & v$nodetypes != "c")
#   edgeintimes <- v$nodetimes[edgesin]
# 
#   # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
#   v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-
#     sample_coaltimes(edgeintimes, tinf.prop, p)
# 
#   # order all edges (transmission, sample, coalescent) by endtime in hostID
#   edgeendtimes <- v$nodetimes[v$nodehosts == hostID]
#   edgeendtypes <- v$nodetypes[v$nodehosts == hostID]
#   nodeorder <- order(edgeendtimes, edgeendtypes)
#   edgesend <- which(v$nodehosts == hostID)[nodeorder]
# 
#   # sample topology of minitree within hostID
#   v$nodeparents[edgesend] <-
#     sample_topology(edgesend,
#                     edgeendtimes[nodeorder],
#                     edgeendtypes[nodeorder],
#                     transnode)
# 
#   copy2pbe1("v", environment())
# }
# 
# rewire_pathCF1_wh_strict <- function() {
#   ### Make input locally available
#   d <- pbe1$d
#   p <- pbe1$p
#   v <- pbe1$v
#   hostID <- pbe1$hostID
# 
#   ### Identify new infector and old infector
#   newinfector <- which(v$infectors == hostID)
#   newinfector <- newinfector[which(v$inftimes[newinfector] == min(v$inftimes[newinfector]))]
#   oldinfector <- v$infectors[hostID]
# 
#   # transmission nodes of host ID and new infector
#   transnode <- 2*d$nsamples - 1 + hostID
#   transnode_ni <- 2*d$nsamples - 1 + newinfector
# 
# 
#   ### Bookkeeping: infection times and transmission tree topology
#   v$nodetimes[c(transnode, transnode_ni)] <- v$nodetimes[c(transnode_ni, transnode)]
#   v$inftimes[c(hostID, newinfector)] <- v$inftimes[c(newinfector, hostID)]
#   v$nodehosts[c(transnode, transnode_ni)] <- c(newinfector, oldinfector)
#   v$infectors[c(hostID, newinfector)] <- c(newinfector, oldinfector)
#   v$nodeparents[transnode_ni] <- v$nodeparents[transnode]
# 
#   ### Move coalescent node from hostID to new infector
#   v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- newinfector
# 
#   ### Change minitrees
#   for(ID in c(hostID, newinfector)) {
#     # times of coalescent events in ID and bottleneck size
#     v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-
#       sample_coaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"],
#                        v$inftimes[ID], p)
# 
#     # order all edges (transmission, sample, coalescent) by endtime in ID
#     edgeendtimes <- v$nodetimes[v$nodehosts == ID]
#     edgeendtypes <- v$nodetypes[v$nodehosts == ID]
#     nodeorder <- order(edgeendtimes, edgeendtypes)
#     edgesend <- which(v$nodehosts == ID)[nodeorder]
# 
#     # sample topology of minitree within ID
#     v$nodeparents[edgesend] <-
#       sample_topology(edgesend,
#                       edgeendtimes[nodeorder],
#                       edgeendtypes[nodeorder],
#                       2 * d$nsamples - 1 + ID)
#   }
#   copy2pbe1("v", environment())
# }
# 
# rewire_pathD_wh_strict <- function() {
#   ### Make input locally available
#   d <- pbe1$d
#   p <- pbe1$p
#   v <- pbe1$v
#   hostID <- pbe1$hostID
#   tinf.prop <- pbe1$tinf.prop
# 
#   ### Identify old index
#   oldindex <- which(v$infectors == 0)
#   transnode_oi <- 2*d$nsamples - 1 + oldindex
#   transnode <- 2*d$nsamples - 1 + hostID
# 
#   ### Bookkeeping: change infection time
#   v$nodetimes[transnode] <- tinf.prop
#   v$inftimes[hostID] <- tinf.prop
# 
#   ### Move coalescent node old infector to hostID
#   # remove from old infector
#   movingcoalnode <- v$nodeparents[transnode]
#   otherchildnode <- setdiff(which(v$nodeparents == movingcoalnode), transnode)
#   v$nodeparents[c(transnode, otherchildnode)] <-
#     c(0L, v$nodeparents[movingcoalnode])
#   # place in hostID
#   v$nodehosts[movingcoalnode] <- hostID
# 
#   ### Bookkeeping: change transmission tree topology
#   v$nodehosts[c(transnode_oi, transnode)] <- c(hostID, 0L)
#   v$infectors[c(oldindex, hostID)] <- c(hostID, 0L)
# 
#   ### Change minitree
#   # edges entering hostID, with endtimes
#   edgesin <- which(v$nodehosts == hostID & v$nodetypes != "c")
#   edgeintimes <- v$nodetimes[edgesin]
# 
#   # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
#   v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-
#     sample_coaltimes(edgeintimes, tinf.prop, p)
# 
#   # order all edges (transmission, sample, coalescent) by endtime in hostID
#   edgeendtimes <- v$nodetimes[v$nodehosts == hostID]
#   edgeendtypes <- v$nodetypes[v$nodehosts == hostID]
#   nodeorder <- order(edgeendtimes, edgeendtypes)
#   edgesend <- which(v$nodehosts == hostID)[nodeorder]
# 
#   # sample topology of minitree within hostID
#   v$nodeparents[edgesend] <-
#     sample_topology(edgesend,
#                     edgeendtimes[nodeorder],
#                     edgeendtypes[nodeorder],
#                     transnode)
# 
#   copy2pbe1("v", environment())
# }
# 
# rewire_pathE_wh_strict <- function() {
#   ### Make input locally available
#   d <- pbe1$d
#   p <- pbe1$p
#   v <- pbe1$v
#   hostID <- pbe1$hostID
#   tinf.prop <- pbe1$tinf.prop
#   infector.proposed.ID <- pbe1$infector.proposed.ID
# 
#   ### Bookkeeping: change infection time
#   transnode <- 2 * d$nsamples - 1 + hostID
#   v$nodetimes[transnode] <- tinf.prop
#   v$inftimes[hostID] <- tinf.prop
# 
#   ### Move coalescent node from old to new infector
#   # remove from old infector
#   movingcoalnode <- v$nodeparents[transnode]
#   otherchildnode <- setdiff(which(v$nodeparents == movingcoalnode), transnode)
#   v$nodeparents[otherchildnode] <- v$nodeparents[movingcoalnode]
#   v$nodehosts[c(transnode, movingcoalnode)] <- -1L
#   # sample time and edge (child node) in new infector
#   v$nodetimes[movingcoalnode] <-
#     sample_singlecoaltime(v$nodetimes[v$nodehosts == infector.proposed.ID & v$nodetypes != "c"],
#                           v$nodetimes[v$nodehosts == infector.proposed.ID & v$nodetypes == "c"],
#                           v$inftimes[hostID],
#                           v$inftimes[infector.proposed.ID],
#                           p)
#   otherchildnode <-
#     sample_singlechildnode(which(v$nodehosts == infector.proposed.ID),
#                            v$nodeparents[v$nodehosts == infector.proposed.ID],
#                            v$nodetimes[v$nodehosts == infector.proposed.ID],
#                            v$nodetimes[movingcoalnode])
#   # place in new infector
#   v$nodeparents[c(movingcoalnode, otherchildnode)] <-
#     c(v$nodeparents[otherchildnode], movingcoalnode)
#   v$nodehosts[c(transnode, movingcoalnode)] <- infector.proposed.ID
# 
#   ### Bookkeeping: change transmission tree topology
#   v$infectors[hostID] <- infector.proposed.ID
# 
#   ### Change minitree
#   # edges entering hostID, with endtimes
#   edgesin <- which(v$nodehosts == hostID & v$nodetypes != "c")
#   edgeintimes <- v$nodetimes[edgesin]
# 
#   # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
#   v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-
#     sample_coaltimes(edgeintimes, tinf.prop, p)
# 
#   # order all edges (transmission, sample, coalescent) by endtime in hostID
#   edgeendtimes <- v$nodetimes[v$nodehosts == hostID]
#   edgeendtypes <- v$nodetypes[v$nodehosts == hostID]
#   nodeorder <- order(edgeendtimes, edgeendtypes)
#   edgesend <- which(v$nodehosts == hostID)[nodeorder]
# 
#   # sample topology of minitree within hostID
#   v$nodeparents[edgesend] <-
#     sample_topology(edgesend,
#                     edgeendtimes[nodeorder],
#                     edgeendtypes[nodeorder],
#                     transnode)
# 
#   copy2pbe1("v", environment())
# }
# 
# rewire_pathCF2_wh_strict <- function() {
#   ### Make input locally available
#   d <- pbe1$d
#   p <- pbe1$p
#   v <- pbe1$v
#   hostID <- pbe1$hostID
#   
#   ### Identify new infector, old infector, and infectees
#   infectees <- which(v$infectors == hostID)
#   newinfector <- infectees[which(v$inftimes[infectees] == min(v$inftimes[infectees]))]
#   oldinfector <- v$infectors[hostID]
#   infectees <- setdiff(infectees, newinfector)
#   infectees_ni <- which(v$infectors == newinfector)
#   
#   # transmission nodes of host ID and new infector
#   transnode <- 2*d$nsamples - 1 + hostID
#   transnode_ni <- 2*d$nsamples - 1 + newinfector
#   
#   # transmission edges entering hostID and new infector, with endtimes
#   edgesin <- setdiff(which(v$nodehosts == hostID & v$nodetypes == "t"),
#                      transnode_ni)
#   edgeintimes <- v$nodetimes[edgesin]
#   edgesin_ni <- which(v$nodehosts == newinfector & v$nodetypes == "t")
#   edgeintimes_ni <- v$nodetimes[edgesin_ni]
#   
#   ### Bookkeeping: infection times, transmission tree topology
#   v$nodetimes[c(transnode, transnode_ni)] <- v$nodetimes[c(transnode_ni, transnode)]
#   v$inftimes[c(hostID, newinfector)] <- v$inftimes[c(newinfector, hostID)]
#   v$nodehosts[c(transnode, edgesin)] <- newinfector
#   v$nodehosts[transnode_ni] <- oldinfector
#   v$nodehosts[edgesin_ni] <- hostID
#   v$infectors[c(hostID, infectees)] <- newinfector
#   v$infectors[newinfector] <- oldinfector
#   v$infectors[infectees_ni] <- hostID
#   
#   ### Swap transmission nodes in phylogenetic tree
#   transnodechild <- which(v$nodeparents == transnode)
#   transnodechild_ni <- which(v$nodeparents == transnode_ni)
#   v$nodeparents[c(transnodechild, transnode, transnodechild_ni, transnode_ni)] <-
#     v$nodeparents[c(transnodechild_ni, transnode_ni, transnodechild, transnode)]
#   
#   ### Move coalescent nodes between hostID and new infector
#   coalnodes <- which(v$nodehosts == hostID & v$nodetypes == "c")
#   coalnodes_ni <- which(v$nodehosts == newinfector & v$nodetypes == "c")
#   # remove all sampling edges from their minitrees
#   sampleedges <- which(v$nodehosts == hostID & (v$nodetypes == "s" | v$nodetypes == "x"))
#   sampleedges_ni <- which(v$nodehosts == newinfector & (v$nodetypes == "s" | v$nodetypes == "x"))
#   matchingcoalnodes <- c()
#   matchingcoalnodes_ni <- c()
#   sampleedge <- sampleedges[3]
#   for(sampleedge in sampleedges) {
#     movingcoalnode <- v$nodeparents[sampleedge]
#     otherchildnode <- setdiff(which(v$nodeparents == movingcoalnode), sampleedge)
#     v$nodeparents[otherchildnode] <- v$nodeparents[movingcoalnode]
#     v$nodeparents[movingcoalnode] <- -1L
#     matchingcoalnodes <- c(matchingcoalnodes, movingcoalnode)
#     coalnodes <- setdiff(coalnodes, movingcoalnode)
#   }
#   for(sampleedge in sampleedges_ni) {
#     movingcoalnode <- v$nodeparents[sampleedge]
#     otherchildnode <- setdiff(which(v$nodeparents == movingcoalnode), sampleedge)
#     if(length(otherchildnode) == 0) {
#       # if new infector did not have infectees, the last sampleedge will end at the transmission node:
#       # in that case, the sampleedge is replaced by a sampleedge from hostID, and one coalescent node
#       # has to move from hostID to the new infector
#       v$nodeparents[sampleedges[1]] <- movingcoalnode
#       edgesin_ni <- sampleedges[1]
#       edgeintimes_ni <- v$nodetimes[edgesin_ni]
#       sampleedges <- sampleedges[-1]
#       v$nodeparents[sampleedge] <- matchingcoalnodes[1]
#       v$nodehosts[matchingcoalnodes[1]] <- newinfector
#       matchingcoalnodes_ni <- c(matchingcoalnodes_ni, matchingcoalnodes[1])
#       matchingcoalnodes <- matchingcoalnodes[-1]
#     } else {
#       v$nodeparents[otherchildnode] <- v$nodeparents[movingcoalnode]
#       v$nodeparents[movingcoalnode] <- -1L
#       matchingcoalnodes_ni <- c(matchingcoalnodes_ni, movingcoalnode)
#       coalnodes_ni <- setdiff(coalnodes_ni, movingcoalnode)
#     }
#   }
#   v$nodehosts[coalnodes] <- newinfector
#   v$nodehosts[coalnodes_ni] <- hostID
#   # then place the sampling edges into their own host
#   if(length(sampleedges) > 0) {
#     for(i in 1:length(sampleedges)) {
#       v$nodetimes[matchingcoalnodes[i]] <-
#         sample_singlecoaltime(c(edgeintimes_ni, head(v$nodetimes[sampleedges], i - 1)),
#                               v$nodetimes[coalnodes_ni],
#                               v$nodetimes[sampleedges[i]],
#                               v$inftimes[hostID], p)
#       otherchildnode <-
#         sample_singlechildnode(c(edgesin_ni, head(sampleedges, i - 1), coalnodes_ni),
#                                v$nodeparents[c(edgesin_ni, head(sampleedges, i - 1), coalnodes_ni)],
#                                v$nodetimes[c(edgesin_ni, head(sampleedges, i - 1), coalnodes_ni)],
#                                v$nodetimes[matchingcoalnodes[i]])
#       v$nodeparents[c(matchingcoalnodes[i], otherchildnode)] <-
#         c(v$nodeparents[otherchildnode], matchingcoalnodes[i])
#       coalnodes_ni <- c(coalnodes_ni, matchingcoalnodes[i])
#     }
#   }
#   for(i in 1:length(sampleedges_ni)) {
#     v$nodetimes[matchingcoalnodes_ni[i]] <-
#       sample_singlecoaltime(c(edgeintimes, v$inftimes[hostID], head(v$nodetimes[sampleedges_ni], i - 1)),
#                             v$nodetimes[coalnodes],
#                             v$nodetimes[sampleedges_ni[i]],
#                             v$inftimes[newinfector], p)
#     otherchildnode <-
#       sample_singlechildnode(c(edgesin, transnode, head(sampleedges_ni, i - 1), coalnodes),
#                              v$nodeparents[c(edgesin, transnode, head(sampleedges_ni, i - 1), coalnodes)],
#                              v$nodetimes[c(edgesin, transnode, head(sampleedges_ni, i - 1), coalnodes)],
#                              v$nodetimes[matchingcoalnodes_ni[i]])
#     v$nodeparents[c(matchingcoalnodes_ni[i], otherchildnode)] <-
#       c(v$nodeparents[otherchildnode], matchingcoalnodes_ni[i])
#     coalnodes <- c(coalnodes, matchingcoalnodes_ni[i])
#   }
#   
#   copy2pbe1("v", environment())
# }
# 
# rewire_pathK_wh_strict <- function() {
#   ### First, dismantle minitree
#   # edges entering hostID, with endtimes
#   edgesin <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes != "c")
#   edgeintimes <- pbe1$v$nodetimes[edgesin]
# 
#   # transmission node of hostID
#   transnode <- 2*pbe1$d$nsamples - 1 + pbe1$hostID
# 
#   # all coalescent nodes in new infector and hostID
#   coalnodes <- which(pbe1$v$nodehosts == pbe1$hostID & pbe1$v$nodetypes == "c")
# 
#   # dismantle topology, move transmission node
#   pbe1$v$nodehosts[coalnodes] <- -1
#   pbe1$v$nodeparents[c(edgesin, coalnodes)] <- -1
# 
#   ### Second, rebuild minitree
#   # times of coalescent events in hostID and bottleneck size, and distribute coalescent nodes over hostID and pre-hostID
#   newcoaltimes <- sample_coaltimes(edgeintimes, pbe1$v$inftimes[pbe1$hostID], pbe1$p)
# 
#   # order all edges (transmission, sample, coalescent) by endtime in hostID
#   nodeorder <- order(c(newcoaltimes, edgeintimes))
#   edgeend <- c(coalnodes, edgesin)[nodeorder]
#   edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]
# 
#   # sample topology of minitree within hostID
#   edgestart <- sample_topology(edgeend,
#                                edgeendtimes,
#                                c(rep("c", length(newcoaltimes)),
#                                  rep("x", length(edgeintimes)))[nodeorder],
#                                transnode)
# 
#   # change minitree in hostID
#   pbe1$v$nodehosts[edgeend] <- pbe1$hostID
#   pbe1$v$nodeparents[edgeend] <- edgestart
#   pbe1$v$nodetimes[edgeend] <- edgeendtimes
# }
# 
# rewire_pullnodes_wh_strict <- function(currentID) {
#   loosenodes <- which(pbe1$v$nodehosts == currentID & pbe1$v$nodeparents == -1)
#   if(length(loosenodes) > 0) {
#     free_cnodes <- which(pbe1$v$nodetypes == "c" & pbe1$v$nodeparents == -1)
#     if(currentID == 0) {
#       pbe1$v$nodeparents[loosenodes] <- 0L
#     } else {
#       # old edges entering currentID, with endtimes
#       edgesendold <- setdiff(which(pbe1$v$nodehosts == currentID & pbe1$v$nodetypes != "c"), loosenodes)
#       if(length(edgesendold) == 0) {
#         parentnode <- 2 * pbe1$d$nsamples - 1 + currentID
#         edgesendold <- currentID
#         pbe1$v$nodeparents[edgesendold] <- parentnode
#         pbe1$v$nodehosts[parentnode] <- pbe1$v$infectors[currentID]
#         pbe1$v$nodetimes[parentnode] <- pbe1$v$inftimes[currentID]
#         loosenodes <- setdiff(loosenodes, currentID)
#       }
#       edgesendoldtimes <- pbe1$v$nodetimes[edgesendold]
# 
#       if(length(loosenodes) > 0) {
#         # coalescentnodes in currentID, with endtimes
#         coalescentnodesold <- which(pbe1$v$nodehosts == currentID & pbe1$v$nodetypes == "c")
#         coalescenttimesold <- pbe1$v$nodetimes[coalescentnodesold]
# 
#         # endtimes of new edges, and new coalescenttimes
#         loosenodetimes <- c()
#         coalescenttimesnew <- c()
#         for(le in 1:length(loosenodes)) {
#           coalescenttimesnew <- c(coalescenttimesnew,
#                                   sample_singlecoaltime(c(edgesendoldtimes, loosenodetimes),
#                                                         c(coalescenttimesold, coalescenttimesnew),
#                                                         pbe1$v$nodetimes[loosenodes[le]], pbe1$v$inftimes[currentID], pbe1$p))
#           loosenodetimes <- c(loosenodetimes, pbe1$v$nodetimes[loosenodes[le]])
#         }
# 
#         # old within-host minitree
#         childnodes <- c(edgesendold, coalescentnodesold)
#         parentnodes <- pbe1$v$nodeparents[childnodes]
#         childnodestimes <- c(edgesendoldtimes, coalescenttimesold)
# 
#         # place new edges into minitree
#         loosenodestoinfector <- c()
#         free_cnodestoinfector <- c()
#         for(le in 1:length(loosenodes)) {
#           newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
#           childnodes <- c(childnodes, loosenodes[le], free_cnodes[le])
#           parentnodes <- c(parentnodes, free_cnodes[le], parentnodes[childnodes == newchildnode])
#           childnodestimes <- c(childnodestimes, loosenodetimes[le], coalescenttimesnew[le])
#           parentnodes[childnodes == newchildnode] <- free_cnodes[le]
#         }
# 
#         # change phybreak object
#         pbe1$v$nodetimes[childnodes] <- childnodestimes
#         pbe1$v$nodehosts[childnodes] <- currentID
#         pbe1$v$nodeparents[childnodes] <- parentnodes
#       }
#     }
#   }
# }
