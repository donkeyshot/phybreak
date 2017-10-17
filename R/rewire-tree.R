# 
# # library(phybreak)
# # set.seed(1)
# # simul <- sim.phybreak(20, wh.model = 2)
# # pbo <- phybreak(simul, wh.model = 2)
# # plotPhylo(pbo)
# # phytree <- get.phylo(pbo)
# # varphybreak <- pbo$v
# # var <- list(
# #   nodetimes = varphybreak$nodetimes[1:39],
# #   nodehosts = varphybreak$nodehosts[1:39],
# #   nodeparents = varphybreak$nodeparents,
# #   nodetypes = varphybreak$nodetypes[1:39],
# #   inftimes = varphybreak$nodetimes[40:59],
# #   infectors = varphybreak$nodehosts[40:59]
# #   )
# # var$nodeparents[var$nodeparents>39] <- (c(0,var$nodeparents)[1+var$nodeparents][var$nodeparents>39])
# # var$nodeparents <- var$nodeparents[1:39]
# # params <- list(wh.model = "infinite", wh.slope = 1, wh.exponent = .3, wh.level = .5, sample.mean = 1)
# # 
# # pbenv <- new.env()
# # list2env(list(v = var, p = params), pbenv)
# 
# # nodechildren <- function(nodeID, var) {
# #   which(var$nodeparents == nodeID)
# # }
# # nodechildhosts <- function(nodeID, var) {
# #   var$nodehosts[nodechildren(nodeID, var)]
# # }
# 
# #all progeny hosts of hostID
# progeny_hosts <- function(infectors, hostID) {
#   which(sapply(1:length(infectors),
#                function(x) hostID %in% c(0, .ptr(infectors, x)[-1])))
# }
# 
# #nodes in hostID or its progeny which have their parent in hostID's ancestors (connected 'edges in bottleneck' of hostID)
# firstnodes <- function(phybreakenv, hostID) {
#   progenyplus <- c(hostID, progeny_hosts(phybreakenv$v$infectors, hostID))
#   which(phybreakenv$v$nodehosts %in% progenyplus & 
#           !(c(0, 0, phybreakenv$v$nodehosts)[2 + phybreakenv$v$nodeparents] %in% progenyplus) &
#     
#           phybreakenv$v$nodeparents != -1)
# }
# 
# #second host in infection chain between first and last
# secondhostinchain <- function(infectors, firsthostID, lasthostID) {
#   tail(setdiff(.ptr(infectors, lasthostID), .ptr(infectors, firsthostID)), 1)
# }
# 
# #nodes in hostID's progeny which have their connected parent in hostID or its ancestors (connected edges entering hostID)
# infectee_nodes <- function(phybreakenv, hostID) {
#   progenyIDs <- progeny_hosts(phybreakenv$v$infectors, hostID)
#   which(phybreakenv$v$nodehosts %in% progenyIDs & 
#           !(c(0, 0, phybreakenv$v$nodehosts)[2 + phybreakenv$v$nodeparents] %in% progenyIDs) &
#           
#           phybreakenv$v$nodeparents != -1 &
#           
#           !(c(0, 0, phybreakenv$v$nodehosts)[2 + phybreakenv$v$nodeparents] == hostID & 
#               c(0, 0, phybreakenv$v$nodeparents)[2 + phybreakenv$v$nodeparents] == -1)
#         )
# }
# 
# ### remove edge in the host of the parent node and all progeny until hostID
# rewire_removeedge <- function(phybreakenv, childnode, hostID = NULL) {
#   if(is.null(hostID)) hostID <- phybreakenv$v$nodehosts[childnode]
#   
#   #find parent node and other child of parent node
#   parentnode <- phybreakenv$v$nodeparents[childnode]
#   if(parentnode != 0) {
#     second_childnode <- setdiff(which(phybreakenv$v$nodeparents == parentnode), childnode)
#     infecteeID <- secondhostinchain(phybreakenv$v$infectors, 
#                                     hostID,
#                                     phybreakenv$v$nodehosts[childnode])
#     
#     #remove parent node and heal the break by linking second child to parent's parent
#     phybreakenv$v$nodeparents[second_childnode] <- phybreakenv$v$nodeparents[parentnode]
#     phybreakenv$v$nodeparents[parentnode] <- -1
#     
#     #time of loose parent edge to latest time possible (transmission time or sampling time)
#     phybreakenv$v$nodetimes[parentnode] <- min(phybreakenv$v$inftimes[infecteeID], phybreakenv$v$nodetimes[childnode])
#     phybreakenv$v$nodehosts[parentnode] <- hostID
#   } 
# }
# 
# rewire_removeinfector <- function(phybreakenv, hostID) {
#   #rootnodes in hostID, or in progeny passing through hostID
#   rootnodes <- firstnodes(phybreakenv, hostID)
# 
#   #remove edges before rootnodes
#   invisible(sapply(rootnodes, rewire_removeedge, phybreakenv = phybreakenv, hostID = phybreakenv$v$infectors[hostID]))
# 
#   #remove infector
#   phybreakenv$v$infectors[hostID] <- -1
#   
#   #remove parentnodes from infector, or ancestors
#   phybreakenv$v$nodehosts[phybreakenv$v$nodeparents[rootnodes]] <- -1
#   phybreakenv$v$nodetimes[phybreakenv$v$nodeparents[rootnodes]] <- phybreakenv$v$inftimes[hostID]
# }
# 
# rewire_stripminitree <- function(phybreakenv, hostID) { 
#   incoming_nodes <- infectee_nodes(phybreakenv, hostID)
#   samplenodes <- which(phybreakenv$v$nodehosts == hostID & phybreakenv$v$nodetypes %in% c("s", "x")) #keep edge ending at "s"
# 
#   invisible(sapply(c(samplenodes, incoming_nodes), rewire_removeedge, phybreakenv = phybreakenv, hostID = hostID))
# }
# 
# rewire_removeinfectiontime <- function(phybreakenv, hostID) {
#   phybreakenv$v$inftimes[hostID] <- -Inf
#   
#   #only remaining edge in hostID, ending at "s", gets length 0
#   phybreakenv$v$nodetimes[phybreakenv$v$nodeparents[hostID]] <- phybreakenv$v$nodetimes[hostID]
# }
# 
# rewire_disconnecthost <- function(phybreakenv, hostID) {
#   rewire_removeinfector(phybreakenv, hostID)
#   rewire_stripminitree(phybreakenv, hostID)
#   rewire_removeinfectiontime(phybreakenv, hostID)
# }
# 
# rewire_assigninfectiontime <- function(phybreakenv, hostID, infectiontime) {
#   #check if minitree has been stripped
#   # if(phybreakenv$v$nodeparents[phybreakenv$v$nodehosts == hostID] != -1) {
#   #   stop("try to assign infection time in connected host")
#   # }
#   phybreakenv$v$inftimes[hostID] <- infectiontime
#   
#   #place parentnode of only edge at infectiontime
#   #phybreakenv$v$nodetimes[phybreakenv$v$nodeparents[hostID]] <- infectiontime
# }
# 
# # rewire_addedge_old <- function(variables, childnode, parameters) {
# #   #determine parentnode and IDs of hosts in infection chain
# #   parentnode <- variables$nodeparents[childnode]
# #   if(parentnode == 0) return()
# #   hostID <- variables$nodehosts[parentnode]
# #   infecteeID <- secondhostinchain(variables$infectors, hostID, variables$nodehosts[childnode])
# # 
# #   #sample coalescent time
# #   newcoalescenttime <- -Inf
# #   while(newcoalescenttime < 0 && hostID != 0) {
# #     hostID <- variables$infectors[infecteeID]
# #     variables$nodehosts[parentnode] <- hostID
# #     
# #     #existing tip nodes and hosts
# #     incomingnodes <- infectee_nodes(variables, hostID)
# #     incominghosts <- unlist(sapply(variables$nodehosts[incomingnodes],
# #                             function(x) if(variables$infectors[x] == hostID) x else Recall(variables$infectors[x])))
# #     samplenodes <- which(variables$nodehosts == hostID & variables$nodeparents != -1 & 
# #                            variables$nodetypes %in% c("s", "x"))
# #     oldtipnodes <- c(incomingnodes, samplenodes)
# #     oldcoalescentnodes <- which(variables$nodehosts == hostID & variables$nodeparents != -1 & variables$nodetypes == "c")
# #     
# #     #existing times
# #     inftimecurrent <- if(hostID != 0) {
# #       variables$inftimes[hostID]
# #     } else {
# #       variables$inftimes[infecteeID] - parameters$sample.mean
# #     }
# #     tiptimes <- c(variables$inftimes[incominghosts], 
# #                        variables$nodetimes[samplenodes]) - inftimecurrent
# #     coalescenttimes <- variables$nodetimes[oldcoalescentnodes] - inftimecurrent
# #     newtiptime <- variables$nodetimes[parentnode] - inftimecurrent
# #     
# #     newcoalescenttime <- sample_singlecoaltime(tiptimes, coalescenttimes, newtiptime, parameters)
# #     
# #     variables$nodetimes[parentnode] <- inftimecurrent
# #     infecteeID <- hostID
# #   }
# #   variables$nodetimes[parentnode] <- variables$nodetimes[parentnode] + newcoalescenttime
# # 
# #   
# #   newchildnode <- sample_singlechildnode(c(oldtipnodes, oldcoalescentnodes), 
# #                                          variables$nodeparents[c(oldtipnodes, oldcoalescentnodes)],
# #                                          variables$nodetimes[c(oldtipnodes, oldcoalescentnodes)],
# #                                          variables$nodetimes[parentnode])
# #   variables$nodeparents[parentnode] <- variables$nodeparents[newchildnode]
# #   variables$nodeparents[newchildnode] <- parentnode
# #   variables$nodehosts[parentnode] <- hostID
# # 
# # }
# 
# rewire_addedge <- function(phybreakenv, parentnode) {
#   #determine parentnode and IDs of hosts in infection chain
#   if(parentnode == 0) return()     #root node
#   hostID <- phybreakenv$v$nodehosts[parentnode]
#   if(hostID == -1) return()     #root node of broken branch
#   if(hostID > 0 && phybreakenv$v$inftimes[hostID] == -Inf) return()     #no known infection time
#   #  infecteeID <- secondhostinchain(phybreakenv$v$infectors, hostID, phybreakenv$v$nodehosts[childnode])
#   
#   #sample coalescent time
#   newcoalescenttime <- -Inf
#   firstround <- TRUE
#   while(newcoalescenttime < 0 && (hostID > 0 || firstround)) {
#     if(!firstround) {
#       #pull edge through hostID and continue in infector
#       phybreakenv$v$nodetimes[parentnode] <- phybreakenv$v$inftimes[hostID]
#       hostID <- phybreakenv$v$infectors[hostID]
#       phybreakenv$v$nodehosts[parentnode] <- hostID
#     }
#     
#     if(hostID == -1 || (hostID > 0 && phybreakenv$v$inftimes[hostID] == -Inf)) {
#       #broken branch or unknown infection time, stop here
#       newcoalescenttime <- 0
#     } else {
#       #existing tip nodes and hosts
#       incomingnodes <- infectee_nodes(phybreakenv, hostID)
#       incominghosts <- unlist(sapply(phybreakenv$v$nodehosts[incomingnodes], secondhostinchain, 
#                               infectors = phybreakenv$v$infectors, firsthostID = hostID))
#       samplenodes <- which(phybreakenv$v$nodehosts == hostID & 
#                              phybreakenv$v$nodetypes %in% c("s", "x") & 
#                              
#                              (c(0, 0, phybreakenv$v$nodeparents)[2 + phybreakenv$v$nodeparents] != -1 |
#                                 c(0, 0, phybreakenv$v$nodehosts)[2 + phybreakenv$v$nodeparents] == -1)
#                            )
#       oldtipnodes <- c(incomingnodes, samplenodes)
#       oldcoalescentnodes <- which(phybreakenv$v$nodehosts == hostID & 
#                                     phybreakenv$v$nodeparents != -1 & phybreakenv$v$nodetypes == "c")
# 
#       #existing times
#       inftimecurrent <- if(hostID != 0) {
#         phybreakenv$v$inftimes[hostID]
#       } else {
#         #one sampling interval back
#         phybreakenv$v$nodetimes[parentnode] - phybreakenv$p$sample.mean
#       }
#       tiptimes <- c(phybreakenv$v$inftimes[incominghosts], 
#                     phybreakenv$v$nodetimes[samplenodes]) - inftimecurrent
#       coalescenttimes <- phybreakenv$v$nodetimes[oldcoalescentnodes] - inftimecurrent
#       newtiptime <- phybreakenv$v$nodetimes[parentnode] - inftimecurrent
#       
#       if(length(tiptimes) > 0) {
#         newcoalescenttime <- sample_singlecoaltime(tiptimes, coalescenttimes, newtiptime, 0, phybreakenv$p)
#       }
#       firstround <- FALSE
#     }
#     
#   }
#   phybreakenv$v$nodetimes[parentnode] <- inftimecurrent + newcoalescenttime
#   
#   if(hostID != -1) {
#     newchildnode <- sample_singlechildnode(c(oldtipnodes, oldcoalescentnodes), 
#                                            phybreakenv$v$nodeparents[c(oldtipnodes, oldcoalescentnodes)],
#                                            c(tiptimes, coalescenttimes),
#                                            newcoalescenttime)
#     phybreakenv$v$nodeparents[parentnode] <- phybreakenv$v$nodeparents[newchildnode]
#     phybreakenv$v$nodeparents[newchildnode] <- parentnode
#   }
# 
# }
# 
# 
# rewire_buildminitree <- function(phybreakenv, hostID) {
#   if(phybreakenv$v$inftimes[hostID] == -Inf) {
#     stop("no infection time known")
#   }
#   
#   #childnodes <- which(phybreakenv$v$nodeparents %in% which(phybreakenv$v$nodeparents == -1 & phybreakenv$v$nodehosts == hostID))
#   parentnodes <- which(phybreakenv$v$nodeparents == -1 & phybreakenv$v$nodehosts == hostID)
#   
#   #invisible(sapply(childnodes, rewire_addedge, phybreakenv = phybreakenv, parameters = parameters))
#   invisible(sapply(parentnodes, rewire_addedge, phybreakenv = phybreakenv))
#   
#   rootnodes <- firstnodes(phybreakenv, hostID)
#   
#   
#   #is there an index case?
#   #no_index <- !any(phybreakenv$v$infectors == 0)
# 
#   #which nodes 
#   # tipnodes <- which(variables$nodeparents == -1 & variables$nodehosts != -1)
#   # incoming_hosts <- with(variables,
#   #                        sapply(nodehosts[tipnodes],
#   #                               function(x) if(x == hostID || infectors[x] == hostID) x else Recall(infectors[x])))
#   # 
#   # tiptimes <- with(variables, 
#   #                  inftimes[incoming_hosts] * (incoming_hosts != hostID) +
#   #                    nodetimes[tipnodes] * (incoming_hosts == hostID))
#   # coalescent_times <- sample_coaltimes(tiptimes - variables$inftimes[hostID], parameters)
#   # coalescent_times <- coalescent_times[coalescent_times >= 0] + variables$inftimes[hostID]
#   # bottleneck <- length(tiptimes) - length(coalescent_times)
#   # alltimes <- c(coalescent_times, tiptimes)
#   # 
#   # coalescent_nodes <- tail(c(c(0)[no_index], which(variables$nodehosts == -1)), -bottleneck)
#   # infector_nodes <- head(c(c(0)[no_index], which(variables$nodehosts == -1)), bottleneck)
#   # variables$nodetimes[coalescent_nodes] <- coalescent_times
#   # 
#   # sortbytime <- order(alltimes)
#   # alltimes <- alltimes[sortbytime]
#   # allnodes <- c(coalescent_nodes, tipnodes)[sortbytime]
#   # alltypes <- c(rep("c", length(coalescent_nodes)), rep("t", length(tipnodes)))[sortbytime]
#   # variables$nodeparents[allnodes] <- sample_topology(allnodes, alltimes, alltypes, infector_nodes)
#   # variables$nodehosts[coalescent_nodes] <- hostID
#   # 
#   # return(variables)
# }
# 
# 
# 
# rewire_assigninfector <- function(phybreakenv, hostID, infectorID) {
#   phybreakenv$v$infectors[hostID] <- infectorID
#   childnodes <- firstnodes(phybreakenv, hostID)
#   parentnodes <- phybreakenv$v$nodeparents[childnodes]
#   
#   if(infectorID == 0 && !(any(parentnodes == 0))) {
#     currentrootnode <- which(phybreakenv$v$nodeparents == 0)
#     phybreakenv$v$nodeparents[currentrootnode] <- parentnodes[1]
#     phybreakenv$v$nodeparents[childnodes[1]] <- 0
#     phybreakenv$v$nodetimes[parentnodes[1]] <- 
#       phybreakenv$v$inftimes[phybreakenv$v$nodehosts[currentrootnode]]
#     phybreakenv$v$nodehosts[parentnodes[1]] <- 
#       phybreakenv$v$infectors[phybreakenv$v$nodehosts[currentrootnode]]
#     parentnodes[1] <- 0
#   }
#   
#   phybreakenv$v$nodehosts[parentnodes] <- infectorID
#   invisible(sapply(parentnodes, rewire_addedge, phybreakenv = phybreakenv))
# }
# 
# rewire_reconnecthost <- function(phybreakenv, hostID, infectorID, infectiontime) {
#   rewire_assigninfectiontime(phybreakenv, hostID, infectiontime)
#   rewire_buildminitree(phybreakenv, hostID)
#   rewire_assigninfector(phybreakenv, hostID, infectorID)
# }
# 
# rewire_changelink <- function(variables, hostID, infectorID, infectiontime) {
#   #rewire_removeinfector
#   #nodetimes: transmission node in infector gets infectiontime
#   #rewire_assigninfector
#   #rewire_assignminitree
# }
# 
# rewire_swaphosts <- function(phybreakenv, hostID, exchange) {
#   #determine infecteeID and change transmission tree
#   infectees <- which(phybreakenv$v$infectors == hostID)
#   infecteeID <- infectees[order(phybreakenv$v$inftimes[infectees])][1]
#   oldinfector <- phybreakenv$v$infectors[hostID]
#   oldfirstinfectiontime <- phybreakenv$v$inftimes[hostID]
#   oldsecondinfectiontime <- phybreakenv$v$inftimes[infecteeID]
# 
#   if(exchange) {
#     ### exchange == TRUE means that infectees and minitrees are exchanged, apart from the sampling nodes
#     #change transmission tree (including exchange of infectees)
#     oldfirstinfectees <- setdiff(infectees, infecteeID)
#     oldsecondinfectees <- which(phybreakenv$v$infectors == infecteeID)
#     phybreakenv$v$inftimes[c(hostID, infecteeID)] <- c(oldsecondinfectiontime, oldfirstinfectiontime)
#     phybreakenv$v$infectors[c(hostID, infecteeID)] <- c(infecteeID, oldinfector)
#     phybreakenv$v$infectors[oldfirstinfectees] <- infecteeID
#     phybreakenv$v$infectors[oldsecondinfectees] <- hostID
#     
#     #remove the samplenodes from the minitrees, so that the minitrees can be exchanged
#     sampleedges <- which(phybreakenv$v$nodehosts %in% c(hostID, infecteeID) &
#                            phybreakenv$v$nodetypes %in% c("s", "x"))
#     invisible(sapply(sampleedges, rewire_removeedge, phybreakenv = phybreakenv))
#     parentnodes <- phybreakenv$v$nodeparents[sampleedges]
#     
#     #exchange the minitrees
#     oldfirstcoalnodes <- setdiff(which(phybreakenv$v$nodehosts == hostID & phybreakenv$v$nodetypes == "c"), parentnodes)
#     oldsecondcoalnodes <- setdiff(which(phybreakenv$v$nodehosts == infecteeID & phybreakenv$v$nodetypes == "c"), parentnodes)
#     phybreakenv$v$nodehosts[oldfirstcoalnodes] <- infecteeID
#     phybreakenv$v$nodehosts[oldsecondcoalnodes] <- hostID
#     
#     #move parent nodes of sample nodes all the way up to the sample nodes, then add them
#     invisible(sapply(parentnodes, rewire_addedge, phybreakenv = phybreakenv))
#     
#   } else {
#     ### exchange == FALSE means that the minitrees are complete resampled
#     #disconnect both hosts
#     rewire_disconnecthost(phybreakenv, hostID)
#     rewire_disconnecthost(phybreakenv, infecteeID)
# 
#     #then reconnect them
#     rewire_reconnecthost(phybreakenv, hostID, infecteeID, oldsecondinfectiontime)
#     rewire_reconnecthost(phybreakenv, infecteeID, oldinfector, oldfirstinfectiontime)
#   }
# }
# 
# 
# 
