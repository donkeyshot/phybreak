# # 
# # library(phybreak)
# # set.seed(1)
# # simul <- sim_phybreak(20, wh.model = 4)
# # pbo <- phybreak(simul, wh.model = 4)
# # plotPhylo(pbo)
# # phytree <- get_phylo(pbo)
# # var <- pbo$v
# # params <- list(obs = 20, wh.model = "exponential", wh.slope = 1, wh.exponent = .3, wh.level = .5, sample.mean = 1, gen.mean = 1)
# # dat <- list(nsamples = 20)
# # 
# # pbenv <- new.env()
# # list2env(list(d = dat, v = var, p = params), pbenv)
# 
# infectoredges <- function(x, hostID, transmissioncluster) {
#   # x is phybreakenvironment
#   which(x$v$nodehosts %in% transmissioncluster & 
#           !(c(0, 0, x$v$nodehosts)[2 + x$v$nodeparents] %in% transmissioncluster))
#   
# }
# 
# infecteeedges <- function(x, hostID, transmissioncluster) {
#   progenyIDs <- setdiff(transmissioncluster, hostID)
#   which(x$v$nodehosts %in% progenyIDs & 
#           !(c(0, 0, x$v$nodehosts)[2 + x$v$nodeparents] %in% progenyIDs))
# }
# 
# ###paths A and E
# rewireEXP_deconrecon <- function(x, hostID, newinfector, newinftime) {
#   # x is phybreakenvironment
#   
#   # transmission edges entering hostID, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs,
#                                       function(y) hostID %in% c(.ptr(x$v$infectors, y))))
#   edgesin <- infecteeedges(x, hostID, transmissioncluster)
#   edgeinhosts <- x$v$nodehosts[edgesin]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgeinhosts %in% infecteeIDs))) {
#     edgeinhosts[tochange] <- x$v$infectors[edgeinhosts][tochange]
#   }
#   edgeintimes <- x$v$inftimes[edgeinhosts]
#   
#   # sample edges in hostID, with endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   
#   # parent nodes of transmission edges and sample edges, consisting of
#   # (1) coalescent nodes in hostID
#   parentnodes <- which(x$v$nodehosts == hostID & x$v$nodetypes == "c")
#   # (2) parents of edges leaving hostID
#   edgesout <- infectoredges(x, hostID, transmissioncluster)
#   for(childnode in edgesout) {
#     parentnode <- x$v$nodeparents[childnode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), childnode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#     }
#     parentnodes <- c(parentnode, parentnodes)
#   }
#   
#   # times of coalescent events in hostID, and bottleneck size
#   newcoaltimes <- newinftime + sample_coaltimes(c(edgeintimes, edgesampletimes) - newinftime, x$p)
#   bottlenecksize <- 1 + sum(newcoaltimes < newinftime)
#   newcoaltimes <- tail(c(0, newcoaltimes), - bottlenecksize)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in hostID, and determine in which host they end
#   nodeorder <- order(c(edgeintimes, edgesampletimes, newcoaltimes))
#   edgeend <- c(edgesin, edgessample, parentnodes[-(1:bottlenecksize)])[nodeorder]
#   edgeendtimes <- c(edgeintimes, edgesampletimes, newcoaltimes)[nodeorder]
#   edgeendhosts <- c(x$v$nodehosts[edgesin], x$v$nodehosts[edgessample], rep(hostID, length(newcoaltimes)))[nodeorder]
#   
#   # sample topology of minitree within hostID
#   edgestart <- sample_topology(edgeend, 
#                                edgeendtimes, 
#                                c(rep("x", length(edgeintimes) + length(edgesampletimes)), 
#                                  rep("c", length(newcoaltimes)))[nodeorder],
#                                parentnodes[1:bottlenecksize])
#   
#   # add edges in the bottleneck going to the infector
#   edgestart <- c(rep(-1, bottlenecksize), edgestart)
#   edgeend <- c(parentnodes[1:bottlenecksize], edgeend)
#   edgeendtimes <- c(rep(newinftime, bottlenecksize), edgeendtimes)
#   edgeendhosts <- c(rep(newinfector, bottlenecksize), edgeendhosts)
#   
#   # in the index case, remove edges ending at node 0 (node 0 is defined as the start)
#   if(any(toremove <- edgeend == 0)) {
#     edgeend <- edgeend[!toremove]
#     edgeendtimes <- edgeendtimes[!toremove]
#     edgeendhosts <- edgeendhosts[!toremove]
#     edgestart <- edgestart[!toremove]
#   }
# 
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeend] <- edgeendhosts
#   x$v$nodeparents[edgeend] <- edgestart
#   x$v$nodetimes[parentnodes[-(1:bottlenecksize)]] <- newcoaltimes
#   x$v$nodetimes[parentnodes[1:bottlenecksize]] <- newinftime
# 
#   # change the transmission tree in the phybreak variables
#   x$v$infectors[hostID] <- newinfector
#   x$v$inftimes[hostID] <- newinftime
#   
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   rewireEXP_connectinfector(x, hostID)
# }
# 
# rewireEXP_connectinfector <- function(x, hostID) {
#   previoushostID <- hostID
#   currenthostID <- x$v$infectors[previoushostID]
#   while(any(x$v$nodeparents == -1) && currenthostID > 0) {
#     # new edges entering currenthostID
#     edgesstartnew <- which(x$v$nodeparents == -1)
#     edgesendnew <- which(x$v$nodeparents %in% edgesstartnew)
#     edgesstartnew <- x$v$nodeparents[edgesendnew]
#     
#     # old edges entering currenthostID, with endtimes
#     transmissioncluster <- which(sapply(1:x$p$obs, 
#                                         function(y) currenthostID %in% c(.ptr(x$v$infectors, y))))
#     edgesendold <- setdiff(infecteeedges(x, currenthostID, transmissioncluster), edgesendnew)
#     edgesendoldhosts <- x$v$nodehosts[edgesendold]
#     infecteeIDs <- which(x$v$infectors == currenthostID)
#     while(any(tochange <- !(edgesendoldhosts %in% infecteeIDs))) {
#       edgesendoldhosts[tochange] <- x$v$infectors[edgesendoldhosts][tochange]
#     }
#     edgesendoldtimes <- x$v$inftimes[edgesendoldhosts]
#     
#     # sample edges in currenthostID, with endtimes
#     edgessample <- which(x$v$nodehosts == currenthostID & x$v$nodetypes %in% c("s", "x"))
#     edgessampletimes <- x$v$nodetimes[edgessample]
#     
#     # coalescentnodes in currenthostID, with endtimes
#     coalescentnodesold <- setdiff(which(x$v$nodehosts == currenthostID & x$v$nodetypes == "c"), edgesstartnew)
#     coalescenttimesold <- x$v$nodetimes[coalescentnodesold]
#     
#     # endtimes of new edges, and new coalescenttimes
#     edgesendnewtimes <- c()
#     coalescenttimesnew <- c()
#     for(le in 1:length(edgesendnew)) {
#       coalescenttimesnew <- c(coalescenttimesnew, 
#                               x$v$inftimes[currenthostID] + 
#                                 sample_singlecoaltime(c(edgesendoldtimes, edgessampletimes, edgesendnewtimes) - x$v$inftimes[currenthostID],
#                                                       c(coalescenttimesold, coalescenttimesnew) - x$v$inftimes[currenthostID],
#                                                       x$v$inftimes[previoushostID] - x$v$inftimes[currenthostID], x$p))
#       edgesendnewtimes <- c(edgesendnewtimes, x$v$inftimes[previoushostID])
#     }
#     
#     # old within-host minitree
#     childnodes <- c(edgesendold, edgessample, coalescentnodesold)
#     parentnodes <- x$v$nodeparents[childnodes]
#     childnodestimes <- c(edgesendoldtimes, edgessampletimes, coalescenttimesold)
#     
#     # place new edges into minitree
#     for(le in 1:length(edgesendnew)) {
#       if(coalescenttimesnew[le] > x$v$inftimes[currenthostID]) {
#         newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
#         childnodes <- c(childnodes, edgesendnew[le], edgesstartnew[le])
#         parentnodes <- c(parentnodes, edgesstartnew[le], parentnodes[childnodes == newchildnode])
#         childnodestimes <- c(childnodestimes, edgesendnewtimes[le], coalescenttimesnew[le])
#         parentnodes[childnodes == newchildnode] <- edgesstartnew[le]
#       }
#     }
#     
#     # change phybreak object
#     x$v$nodetimes[edgesstartnew] <- coalescenttimesnew
#     x$v$nodeparents[childnodes] <- parentnodes
#     x$v$nodehosts[x$v$nodeparents == -1] <- x$v$infectors[currenthostID]
#     x$v$nodetimes[x$v$nodeparents == -1] <- x$v$inftimes[currenthostID]
#     
#     # go to infector
#     previoushostID <- currenthostID
#     currenthostID <- x$v$infectors[previoushostID]
#   }
#   
#   if(any(x$v$nodeparents == -1) && currenthostID == 0) rewireEXP_connectindex(x, previoushostID)
# }
# 
# rewireEXP_connectindex <- function(x, hostID) {
#   # new edges leaving index
#   edgesstartnew <- which(x$v$nodeparents == -1 & x$v$nodehosts == 0)
#   edgesendnew <- which(x$v$nodeparents %in% edgesstartnew)
#   edgesstartnew <- x$v$nodeparents[edgesendnew]
#   
#   # old edges leaving index, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs, 
#                                       function(y) hostID %in% c(.ptr(x$v$infectors, y))))
#   edgesendold <- setdiff(infecteeedges(x, 0, transmissioncluster), edgesendnew)
#   edgesendoldhosts <- rep(hostID, length(edgesendold))
#   edgesendoldtimes <- x$v$inftimes[edgesendoldhosts]
#   
#   # coalescentnodes before index, with endtimes
#   coalescentnodesold <- setdiff(which(x$v$nodehosts == 0 & x$v$nodetypes == "c"), edgesstartnew)
#   coalescenttimesold <- x$v$nodetimes[coalescentnodesold]
#   
#   # endtimes of new edges, and new coalescenttimes
#   edgesendnewtimes <- c()
#   coalescenttimesnew <- c()
#   if(length(edgesendnew) > 0) {
#     for(le in 1:length(edgesendnew)) {
#       coalescenttimesnew <- c(coalescenttimesnew, x$v$inftimes[hostID] - x$p$gen.mean + 
#                                 sample_singlecoaltime(c(edgesendoldtimes, edgesendnewtimes) - x$v$inftimes[hostID] + x$p$gen.mean,
#                                                       c(coalescenttimesold, coalescenttimesnew) - x$v$inftimes[hostID] + x$p$gen.mean,
#                                                       x$p$gen.mean, x$p))
#       edgesendnewtimes <- c(edgesendnewtimes, x$v$inftimes[hostID])
#     }
#   }
#   
#   # old within-host minitree
#   childnodes <- c(edgesendold, coalescentnodesold)
#   parentnodes <- x$v$nodeparents[childnodes]
#   childnodestimes <- c(edgesendoldtimes, coalescenttimesold)
#   
#   # place new edges into minitree
#   if(length(edgesendnew) > 0) {
#     for(le in 1:length(edgesendnew)) {
#       newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
#       childnodes <- c(childnodes, edgesendnew[le], edgesstartnew[le])
#       parentnodes <- c(parentnodes, edgesstartnew[le], parentnodes[childnodes == newchildnode])
#       childnodestimes <- c(childnodestimes, edgesendnewtimes[le], coalescenttimesnew[le])
#       parentnodes[childnodes == newchildnode] <- edgesstartnew[le]
#     }
#   }
#   
#   # change phybreak object
#   x$v$nodetimes[edgesstartnew] <- coalescenttimesnew
#   x$v$nodeparents[childnodes] <- parentnodes
# }
# 
# # x <- pbenv
# # hostID <- 2
# # newinftime <- -1
# 
# ### NEW INDEX MUST FIRST BE REMOVED FROM PHYLOTREE!!
# 
# ###path D
# rewireEXP_becomeindex <- function(x, hostID, newinftime) {
#   # x is phybreakenvironment
#   oldindexID <- which(x$v$infectors == 0)
#   
#   # transmission edges entering hostID, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs,
#                                       function(y) hostID %in% c(.ptr(x$v$infectors, y))))
#   edgesin <- infecteeedges(x, hostID, transmissioncluster)
#   edgeinhosts <- x$v$nodehosts[edgesin]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgeinhosts %in% infecteeIDs))) {
#     edgeinhosts[tochange] <- x$v$infectors[edgeinhosts][tochange]
#   }
#   edgeintimes <- x$v$inftimes[edgeinhosts]
#   
#   # sample edges in hostID, with endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   
#   # parent nodes of transmission edges and sample edges, consisting of
#   # (1) coalescent nodes in hostID
#   parentnodes <- which(x$v$nodehosts == hostID & x$v$nodetypes == "c")
#   # (2) parents of edges leaving hostID
#   edgesout <- infectoredges(x, hostID, transmissioncluster)
#   for(childnode in edgesout) {
#     parentnode <- x$v$nodeparents[childnode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), childnode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#       x$v$nodehosts[parentnode] <- hostID
#     }
#     parentnodes <- c(parentnode, parentnodes)
#   }
#   
#   # edges leaving oldindexID
#   transmissionclusterindex <- setdiff(1:x$p$obs, transmissioncluster)
#   edgesoutindex <- infectoredges(x, oldindexID, transmissionclusterindex)
#   parentnodesindex <- c(0, which(x$v$nodehosts == 0))
#   
#   # add old index to edges entering hostID
#   edgesin <- c(edgesin, edgesoutindex)
#   edgeinhosts <- c(edgeinhosts, oldindexID)
#   edgeintimes <- c(edgeintimes, x$v$inftimes[oldindexID])
#   
#   # (3) parent nodes of old index
#   parentnodes <- c(parentnodesindex, parentnodes)
#     
#   # # dismantle minitree before oldindexID
#   # x$v$nodehosts[edgesoutindex] <- hostID
#   # x$v$nodeparents[edgesoutindex] <- parentnodesindex
#   # x$v$nodeparents[parentnodesindex] <- -1
#   # x$v$nodetimes[parentnodesindex] <- x$v$inftimes[oldindexID]
#   
#   # change transmission tree
#   x$v$infectors[hostID] <- 0
#   x$v$infectors[oldindexID] <- hostID
#   x$v$inftimes[hostID] <- newinftime
#   
# 
#   
#   # times of coalescent events in hostID, and bottleneck size
#   newcoaltimes <- newinftime + sample_coaltimes(c(edgeintimes, edgesampletimes) - newinftime, x$p)
#   bottlenecksize <- 1 + sum(newcoaltimes < newinftime)
#   newcoaltimes <- tail(c(0, newcoaltimes), - bottlenecksize)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in hostID, and determine in which host they end
#   nodeorder <- order(c(edgeintimes, edgesampletimes, newcoaltimes))
#   edgeend <- c(edgesin, edgessample, parentnodes[-(1:bottlenecksize)])[nodeorder]
#   edgeendtimes <- c(edgeintimes, edgesampletimes, newcoaltimes)[nodeorder]
#   edgeendhosts <- c(x$v$nodehosts[edgesin], x$v$nodehosts[edgessample], rep(hostID, length(newcoaltimes)))[nodeorder]
#   
#   # sample topology of minitree within hostID
#   edgestart <- sample_topology(edgeend, 
#                                edgeendtimes, 
#                                c(rep("x", length(edgeintimes) + length(edgesampletimes)), 
#                                  rep("c", length(newcoaltimes)))[nodeorder],
#                                parentnodes[1:bottlenecksize])
#   
#   # add edges in the bottleneck going before the index
#   edgestart <- c(rep(-1, bottlenecksize), edgestart)
#   edgeend <- c(parentnodes[1:bottlenecksize], edgeend)
#   edgeendtimes <- c(rep(newinftime, bottlenecksize), edgeendtimes)
#   edgeendhosts <- c(rep(0, bottlenecksize), edgeendhosts)
#   
#   # as it is now the index case, remove edges ending at node 0 (node 0 is defined as the start)
#   toremove <- edgeend == 0
#   edgeend <- edgeend[!toremove]
#   edgeendtimes <- edgeendtimes[!toremove]
#   edgeendhosts <- edgeendhosts[!toremove]
#   edgestart <- edgestart[!toremove]
#   
#   
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeend] <- edgeendhosts
#   x$v$nodeparents[edgeend] <- edgestart
#   x$v$nodetimes[parentnodes[-(1:bottlenecksize)]] <- newcoaltimes
#   x$v$nodetimes[parentnodes[1:bottlenecksize]] <- newinftime
# 
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   # rewireEXP_connectinfector(x, hostID)
#   rewireEXP_connectindex(x, hostID)
# }
# 
# ###path B  
# rewireEXP_resignindex <- function(x, hostID, newinfector, newinftime) {
#   # x is phybreakenvironment
# 
#   # identify new index
#   infecteeIDs <- which(x$v$infectors == hostID)
#   newindex <- infecteeIDs[which(x$v$inftimes[infecteeIDs] == min(x$v$inftimes[infecteeIDs]))]
# 
#   # transmission edges entering hostID, with endtimes
#   transmissioncluster <- 1:x$p$obs
#   edgesin <- infecteeedges(x, hostID, transmissioncluster)
#   edgeinhosts <- x$v$nodehosts[edgesin]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgeinhosts %in% infecteeIDs))) {
#     edgeinhosts[tochange] <- x$v$infectors[edgeinhosts][tochange]
#   }
#   edgeintimes <- x$v$inftimes[edgeinhosts]
#   
#   # identify edges leaving new index and remove them from old index
#   toremove <- edgeinhosts == newindex
#   edgesoutnewindex <- edgesin[toremove]
#   edgesin <- edgesin[!toremove]
#   edgeinhosts <- edgeinhosts[!toremove]
#   edgeintimes <- edgeintimes[!toremove]
#   
#   # sample edges in hostID, with endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   
#   # parent nodes of transmission edges and sample edges, consisting of
#   # (1) coalescent nodes in hostID
#   parentnodes <- which(x$v$nodehosts == hostID & x$v$nodetypes == "c")
#   # (2) coalescent nodes before hostID
#   parentnodes <- c(0, which(x$v$nodehosts == 0), parentnodes)
# 
#   
#   # edges leaving newindex
#   # transmissionclusternewindex <- which(sapply(1:x$p$obs,
#   #                                             function(y) newindex %in% c(.ptr(x$v$infectors, y))))
#   # edgesoutnewindex <- infectoredges(x, newindex, transmissionclusternewindex)
#   
#   # distribute parent nodes over hostID and edges leaving new index
#   bottlenecksizenewindex <- length(edgesoutnewindex)
#   parentnodesnewindex <- parentnodes[1:bottlenecksizenewindex]
#   parentnodes <- parentnodes[-(1:bottlenecksizenewindex)]
#   
#   # dismantle minitree in hostID and move parent nodes to new index
#   x$v$nodehosts[parentnodesnewindex] <- 0
#   x$v$nodeparents[edgesoutnewindex] <- parentnodesnewindex
#   x$v$nodeparents[parentnodesnewindex] <- -1
#   x$v$nodetimes[parentnodesnewindex] <- x$v$inftimes[newindex]
#   x$v$nodehosts[parentnodes] <- hostID
#   x$v$nodeparents[c(edgesin, edgessample)] <- parentnodes
#   x$v$nodeparents[parentnodes] <- -1
#   x$v$nodetimes[parentnodes] <- newinftime
#   
#   # build minitree before new index
#   rewireEXP_connectindex(x, newindex)
# 
#   # times of coalescent events in hostID, and bottleneck size
#   newcoaltimes <- newinftime + sample_coaltimes(c(edgeintimes, edgesampletimes) - newinftime, x$p)
#   bottlenecksize <- 1 + sum(newcoaltimes < newinftime)
#   newcoaltimes <- tail(c(0, newcoaltimes), - bottlenecksize)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in hostID, and determine in which host they end
#   nodeorder <- order(c(edgeintimes, edgesampletimes, newcoaltimes))
#   edgeend <- c(edgesin, edgessample, parentnodes[-(1:bottlenecksize)])[nodeorder]
#   edgeendtimes <- c(edgeintimes, edgesampletimes, newcoaltimes)[nodeorder]
#   edgeendhosts <- c(x$v$nodehosts[edgesin], x$v$nodehosts[edgessample], rep(hostID, length(newcoaltimes)))[nodeorder]
#   
#   # sample topology of minitree within hostID
#   edgestart <- sample_topology(edgeend, 
#                                edgeendtimes, 
#                                c(rep("x", length(edgeintimes) + length(edgesampletimes)), 
#                                  rep("c", length(newcoaltimes)))[nodeorder],
#                                parentnodes[1:bottlenecksize])
#   
#   # add edges in the bottleneck going to the infector
#   edgestart <- c(rep(-1, bottlenecksize), edgestart)
#   edgeend <- c(parentnodes[1:bottlenecksize], edgeend)
#   edgeendtimes <- c(rep(newinftime, bottlenecksize), edgeendtimes)
#   edgeendhosts <- c(rep(newinfector, bottlenecksize), edgeendhosts)
#   
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeend] <- edgeendhosts
#   x$v$nodeparents[edgeend] <- edgestart
#   x$v$nodetimes[parentnodes[-(1:bottlenecksize)]] <- newcoaltimes
#   x$v$nodetimes[parentnodes[1:bottlenecksize]] <- newinftime
#   
#   # change transmission tree
#   x$v$infectors[hostID] <- newindex
#   x$v$infectors[newindex] <- 0
#   x$v$inftimes[hostID] <- newinftime
# 
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   rewireEXP_connectinfector(x, hostID)
# }
#   
# ###path C1
# rewireEXP_swapincompleteindex <- function(x, hostID) {
#   # x is phybreakenvironment
#   
#   # identify new index and new infection times
#   infecteeIDs <- which(x$v$infectors == hostID)
#   newindex <- infecteeIDs[which(x$v$inftimes[infecteeIDs] == min(x$v$inftimes[infecteeIDs]))]
#   newinftime_hostID <- x$v$inftimes[newindex]
#   newinftime_newindex <- x$v$inftimes[hostID]
#   
#   # transmission edges entering hostID, with endtimes
#   transmissioncluster <- 1:x$p$obs
#   edgesin <- infecteeedges(x, hostID, transmissioncluster)
#   edgeinhosts <- x$v$nodehosts[edgesin]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgeinhosts %in% infecteeIDs))) {
#     edgeinhosts[tochange] <- x$v$infectors[edgeinhosts][tochange]
#   }
#   edgeintimes <- x$v$inftimes[edgeinhosts]
#   
#   # identify edges leaving new index and remove them from old index
#   toremove <- edgeinhosts == newindex
#   edgesin <- edgesin[!toremove]
#   edgeinhosts <- edgeinhosts[!toremove]
#   edgeintimes <- edgeintimes[!toremove]
#   
#   # sample edges in hostID, with endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   
#   # transmission edges entering new index, with endtimes
#   transmissionclusternewindex <- which(sapply(1:x$p$obs,
#                                               function(y) newindex %in% c(.ptr(x$v$infectors, y))))
#   edgesinnewindex <- infecteeedges(x, newindex, transmissionclusternewindex)
#   edgeinhostsnewindex <- x$v$nodehosts[edgesinnewindex]
#   infecteeIDs <- which(x$v$infectors == newindex)
#   while(any(tochange <- !(edgeinhostsnewindex %in% infecteeIDs))) {
#     edgeinhostsnewindex[tochange] <- x$v$infectors[edgeinhostsnewindex][tochange]
#   }
#   edgeintimesnewindex <- x$v$inftimes[edgeinhostsnewindex]
#   
#   # sample edges in new index, with endtimes
#   edgessamplenewindex <- which(x$v$nodehosts == newindex & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimesnewindex <- x$v$nodetimes[edgessamplenewindex]
#   
#   
#   # parent nodes of transmission edges and sample edges in hostID and new index, consisting of
#   # (1) coalescent nodes in hostID and new index
#   parentnodes <- which(x$v$nodehosts %in% c(hostID, newindex) & x$v$nodetypes == "c")
#   # (2) coalescent nodes before hostID
#   parentnodes <- c(0, which(x$v$nodehosts == 0), parentnodes)
#   
#   # distribute parent nodes over hostID and edges leaving new index
#   edgenumbernewindex <- length(c(edgesinnewindex, edgessamplenewindex))
#   parentnodesnewindex <- parentnodes[1:edgenumbernewindex]
#   parentnodes <- parentnodes[-(1:edgenumbernewindex)]
#   
#   # times of coalescent events in hostID, and bottleneck size
#   newcoaltimes <- newinftime_hostID + sample_coaltimes(c(edgeintimes, edgesampletimes) - newinftime_hostID, x$p)
#   bottlenecksize <- 1 + sum(newcoaltimes < newinftime_hostID)
#   newcoaltimes <- tail(c(0, newcoaltimes), - bottlenecksize)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in hostID, and determine in which host they end
#   nodeorder <- order(c(edgeintimes, edgesampletimes, newcoaltimes))
#   edgeend <- c(edgesin, edgessample, parentnodes[-(1:bottlenecksize)])[nodeorder]
#   edgeendtimes <- c(edgeintimes, edgesampletimes, newcoaltimes)[nodeorder]
#   edgeendhosts <- c(x$v$nodehosts[edgesin], x$v$nodehosts[edgessample], rep(hostID, length(newcoaltimes)))[nodeorder]
#   
#   # sample topology of minitree within hostID
#   edgestart <- sample_topology(edgeend, 
#                                edgeendtimes, 
#                                c(rep("x", length(edgeintimes) + length(edgesampletimes)), 
#                                  rep("c", length(newcoaltimes)))[nodeorder],
#                                parentnodes[1:bottlenecksize])
#   
#   # add edges in the bottleneck going to the infector
#   edgestart <- c(rep(-1, bottlenecksize), edgestart)
#   edgeend <- c(parentnodes[1:bottlenecksize], edgeend)
#   edgeendtimes <- c(rep(newinftime_hostID, bottlenecksize), edgeendtimes)
#   edgeendhosts <- c(rep(newindex, bottlenecksize), edgeendhosts)
#   
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeend] <- edgeendhosts
#   x$v$nodeparents[edgeend] <- edgestart
#   x$v$nodetimes[parentnodes[-(1:bottlenecksize)]] <- newcoaltimes
#   x$v$nodetimes[parentnodes[1:bottlenecksize]] <- newinftime_hostID
#   
#   # add bottleneck edges to edges going into new index
#   parentnodesnewindex <- c(parentnodesnewindex, parentnodes[1:bottlenecksize])
#   edgesinnewindex <- c(edgesinnewindex, edgeend[edgestart %in% parentnodes[1:bottlenecksize]])
#   edgeinhostsnewindex <- c(edgeinhostsnewindex, rep(hostID, bottlenecksize))
#   edgeintimesnewindex <- c(edgeintimesnewindex, rep(newinftime_hostID, bottlenecksize))
#   
#   # times of coalescent events in new index, and bottleneck size
#   newcoaltimesnewindex <- newinftime_newindex + 
#     sample_coaltimes(c(edgeintimesnewindex, edgesampletimesnewindex) - newinftime_newindex, x$p)
#   bottlenecksizenewindex <- 1 + sum(newcoaltimesnewindex < newinftime_newindex)
#   newcoaltimesnewindex <- tail(c(0, newcoaltimesnewindex), - bottlenecksizenewindex)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in hostID, and determine in which host they end
#   nodeordernewindex <- order(c(edgeintimesnewindex, edgesampletimesnewindex, newcoaltimesnewindex))
#   edgeendnewindex <- c(edgesinnewindex, edgessamplenewindex, parentnodesnewindex[-(1:bottlenecksizenewindex)])[nodeordernewindex]
#   edgeendtimesnewindex <- c(edgeintimesnewindex, edgesampletimesnewindex, newcoaltimesnewindex)[nodeordernewindex]
#   edgeendhostsnewindex <- c(x$v$nodehosts[edgesinnewindex], x$v$nodehosts[edgessamplenewindex], 
#                     rep(newindex, length(newcoaltimesnewindex)))[nodeordernewindex]
#   
#   # sample topology of minitree within hostID
#   edgestartnewindex <- sample_topology(edgeendnewindex, 
#                                        edgeendtimesnewindex, 
#                                c(rep("x", length(edgeintimesnewindex) + length(edgesampletimesnewindex)), 
#                                  rep("c", length(newcoaltimesnewindex)))[nodeordernewindex],
#                                parentnodesnewindex[1:bottlenecksizenewindex])
#   
#   # add edges in the bottleneck going to the infector
#   edgestartnewindex <- c(rep(-1, bottlenecksizenewindex), edgestartnewindex)
#   edgeendnewindex <- c(parentnodesnewindex[1:bottlenecksizenewindex], edgeendnewindex)
#   edgeendtimesnewindex <- c(rep(newinftime_newindex, bottlenecksizenewindex), edgeendtimesnewindex)
#   edgeendhostsnewindex <- c(rep(0, bottlenecksizenewindex), edgeendhostsnewindex)
#   
#   # as it is now the index case, remove edges ending at node 0 (node 0 is defined as the start)
#   toremove <- edgeendnewindex == 0
#   edgeendnewindex <- edgeendnewindex[!toremove]
#   edgeendtimesnewindex <- edgeendtimesnewindex[!toremove]
#   edgeendhostsnewindex <- edgeendhostsnewindex[!toremove]
#   edgestartnewindex <- edgestartnewindex[!toremove]
#   
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeendnewindex] <- edgeendhostsnewindex
#   x$v$nodeparents[edgeendnewindex] <- edgestartnewindex
#   x$v$nodetimes[parentnodesnewindex[-(1:bottlenecksizenewindex)]] <- newcoaltimesnewindex
#   x$v$nodetimes[parentnodesnewindex[1:bottlenecksizenewindex]] <- newinftime_newindex
#   
# 
#   # change transmission tree
#   x$v$infectors[hostID] <- newindex
#   x$v$infectors[newindex] <- 0
#   x$v$inftimes[hostID] <- newinftime_hostID
#   x$v$inftimes[newindex] <- newinftime_newindex
#   
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   rewireEXP_connectindex(x, newindex)
# }
# 
# ###path F1
# rewireEXP_swapincomplete <- function(x, hostID) {
#   # x is phybreakenvironment
#   
#   # identify new infector and new infection times
#   infectorID <- x$v$infectors[hostID]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   newinfector <- infecteeIDs[which(x$v$inftimes[infecteeIDs] == min(x$v$inftimes[infecteeIDs]))]
#   newinftime_hostID <- x$v$inftimes[newinfector]
#   newinftime_newinfector <- x$v$inftimes[hostID]
#   
#   # transmission edges entering hostID, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs,
#                                       function(y) hostID %in% c(.ptr(x$v$infectors, y))))
#   edgesin <- infecteeedges(x, hostID, transmissioncluster)
#   edgeinhosts <- x$v$nodehosts[edgesin]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgeinhosts %in% infecteeIDs))) {
#     edgeinhosts[tochange] <- x$v$infectors[edgeinhosts][tochange]
#   }
#   edgeintimes <- x$v$inftimes[edgeinhosts]
#   
#   # identify edges leaving new infector and remove them from old infector
#   toremove <- edgeinhosts == newinfector
#   edgesin <- edgesin[!toremove]
#   edgeinhosts <- edgeinhosts[!toremove]
#   edgeintimes <- edgeintimes[!toremove]
#   
#   # sample edges in hostID, with endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   
#   # transmission edges entering new infector, with endtimes
#   transmissionclusternewinfector <- which(sapply(1:x$p$obs,
#                                               function(y) newinfector %in% c(.ptr(x$v$infectors, y))))
#   edgesinnewinfector <- infecteeedges(x, newinfector, transmissionclusternewinfector)
#   edgeinhostsnewinfector <- x$v$nodehosts[edgesinnewinfector]
#   infecteeIDs <- which(x$v$infectors == newinfector)
#   while(any(tochange <- !(edgeinhostsnewinfector %in% infecteeIDs))) {
#     edgeinhostsnewinfector[tochange] <- x$v$infectors[edgeinhostsnewinfector][tochange]
#   }
#   edgeintimesnewinfector <- x$v$inftimes[edgeinhostsnewinfector]
#   
#   # sample edges in new infector, with endtimes
#   edgessamplenewinfector <- which(x$v$nodehosts == newinfector & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimesnewinfector <- x$v$nodetimes[edgessamplenewinfector]
#   
#   # parent nodes of transmission edges and sample edges in hostID and new infector, consisting of
#   # (1) coalescent nodes in hostID and new infector
#   parentnodes <- which(x$v$nodehosts %in% c(hostID, newinfector) & x$v$nodetypes == "c")
#   # (2) coalescent nodes at start of edges leaving hostID
#   edgesout <- infectoredges(x, hostID, transmissioncluster)
#   for(childnode in edgesout) {
#     parentnode <- x$v$nodeparents[childnode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), childnode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#     }
#     parentnodes <- c(parentnode, parentnodes)
#   }
#   
#   # distribute parent nodes over hostID and edges leaving new infector
#   edgenumbernewinfector <- length(c(edgesinnewinfector, edgessamplenewinfector))
#   parentnodesnewinfector <- parentnodes[1:edgenumbernewinfector]
#   parentnodes <- parentnodes[-(1:edgenumbernewinfector)]
#   
#   # times of coalescent events in hostID, and bottleneck size
#   newcoaltimes <- newinftime_hostID + sample_coaltimes(c(edgeintimes, edgesampletimes) - newinftime_hostID, x$p)
#   bottlenecksize <- 1 + sum(newcoaltimes < newinftime_hostID)
#   newcoaltimes <- tail(c(0, newcoaltimes), - bottlenecksize)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in hostID, and determine in which host they end
#   nodeorder <- order(c(edgeintimes, edgesampletimes, newcoaltimes))
#   edgeend <- c(edgesin, edgessample, parentnodes[-(1:bottlenecksize)])[nodeorder]
#   edgeendtimes <- c(edgeintimes, edgesampletimes, newcoaltimes)[nodeorder]
#   edgeendhosts <- c(x$v$nodehosts[edgesin], x$v$nodehosts[edgessample], rep(hostID, length(newcoaltimes)))[nodeorder]
#   
#   # sample topology of minitree within hostID
#   edgestart <- sample_topology(edgeend, 
#                                edgeendtimes, 
#                                c(rep("x", length(edgeintimes) + length(edgesampletimes)), 
#                                  rep("c", length(newcoaltimes)))[nodeorder],
#                                parentnodes[1:bottlenecksize])
#   
#   # add edges in the bottleneck going to the infector
#   edgestart <- c(rep(-1, bottlenecksize), edgestart)
#   edgeend <- c(parentnodes[1:bottlenecksize], edgeend)
#   edgeendtimes <- c(rep(newinftime_hostID, bottlenecksize), edgeendtimes)
#   edgeendhosts <- c(rep(newinfector, bottlenecksize), edgeendhosts)
#   
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeend] <- edgeendhosts
#   x$v$nodeparents[edgeend] <- edgestart
#   x$v$nodetimes[parentnodes[-(1:bottlenecksize)]] <- newcoaltimes
#   x$v$nodetimes[parentnodes[1:bottlenecksize]] <- newinftime_hostID
#   
#   # add bottleneck edges to edges going into new infector
#   parentnodesnewinfector <- c(parentnodesnewinfector, parentnodes[1:bottlenecksize])
#   edgesinnewinfector <- c(edgesinnewinfector, edgeend[edgestart %in% parentnodes[1:bottlenecksize]])
#   edgeinhostsnewinfector <- c(edgeinhostsnewinfector, rep(hostID, bottlenecksize))
#   edgeintimesnewinfector <- c(edgeintimesnewinfector, rep(newinftime_hostID, bottlenecksize))
#   
#   # times of coalescent events in new infector, and bottleneck size
#   newcoaltimesnewinfector <- newinftime_newinfector + 
#     sample_coaltimes(c(edgeintimesnewinfector, edgesampletimesnewinfector) - newinftime_newinfector, x$p)
#   bottlenecksizenewinfector <- 1 + sum(newcoaltimesnewinfector < newinftime_newinfector)
#   newcoaltimesnewinfector <- tail(c(0, newcoaltimesnewinfector), - bottlenecksizenewinfector)
#   
#   # order all edges (transmission, sample, coalescent) by endtime in new infector, and determine in which host they end
#   nodeordernewinfector <- order(c(edgeintimesnewinfector, edgesampletimesnewinfector, newcoaltimesnewinfector))
#   edgeendnewinfector <- c(edgesinnewinfector, edgessamplenewinfector, parentnodesnewinfector[-(1:bottlenecksizenewinfector)])[nodeordernewinfector]
#   edgeendtimesnewinfector <- c(edgeintimesnewinfector, edgesampletimesnewinfector, newcoaltimesnewinfector)[nodeordernewinfector]
#   edgeendhostsnewinfector <- c(x$v$nodehosts[edgesinnewinfector], x$v$nodehosts[edgessamplenewinfector], 
#                             rep(newinfector, length(newcoaltimesnewinfector)))[nodeordernewinfector]
#   
#   # sample topology of minitree within new infector
#   edgestartnewinfector <- sample_topology(edgeendnewinfector, 
#                                        edgeendtimesnewinfector, 
#                                        c(rep("x", length(edgeintimesnewinfector) + length(edgesampletimesnewinfector)), 
#                                          rep("c", length(newcoaltimesnewinfector)))[nodeordernewinfector],
#                                        parentnodesnewinfector[1:bottlenecksizenewinfector])
#   
#   # add edges in the bottleneck going to the infector's infector
#   edgestartnewinfector <- c(rep(-1, bottlenecksizenewinfector), edgestartnewinfector)
#   edgeendnewinfector <- c(parentnodesnewinfector[1:bottlenecksizenewinfector], edgeendnewinfector)
#   edgeendtimesnewinfector <- c(rep(newinftime_newinfector, bottlenecksizenewinfector), edgeendtimesnewinfector)
#   edgeendhostsnewinfector <- c(rep(infectorID, bottlenecksizenewinfector), edgeendhostsnewinfector)
#   
#   # remove edges ending at node 0 (node 0 is defined as the start)
#   if(any(toremove <- edgeendnewinfector == 0)) {
#     edgeendnewinfector <- edgeendnewinfector[!toremove]
#     edgeendtimesnewinfector <- edgeendtimesnewinfector[!toremove]
#     edgeendhostsnewinfector <- edgeendhostsnewinfector[!toremove]
#     edgestartnewinfector <- edgestartnewinfector[!toremove]
#   }
#   
#   # change the minitree in the phybreak variables
#   x$v$nodehosts[edgeendnewinfector] <- edgeendhostsnewinfector
#   x$v$nodeparents[edgeendnewinfector] <- edgestartnewinfector
#   x$v$nodetimes[parentnodesnewinfector[-(1:bottlenecksizenewinfector)]] <- newcoaltimesnewinfector
#   x$v$nodetimes[parentnodesnewinfector[1:bottlenecksizenewinfector]] <- newinftime_newinfector
#   
#   
#   # change transmission tree
#   x$v$infectors[hostID] <- newinfector
#   x$v$infectors[newinfector] <- infectorID
#   x$v$inftimes[hostID] <- newinftime_hostID
#   x$v$inftimes[newinfector] <- newinftime_newinfector
#   
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   rewireEXP_connectinfector(x, newinfector)
# }
# 
# 
# ###path C2
# rewireEXP_swapcompleteindex <- function(x, hostID) {
#   # x is phybreakenvironment
#   
#   # identify new index and new infection times
#   infecteeIDs <- which(x$v$infectors == hostID)
#   newindex <- infecteeIDs[which(x$v$inftimes[infecteeIDs] == min(x$v$inftimes[infecteeIDs]))]
#   newinftime_hostID <- x$v$inftimes[newindex]
#   newinftime_newindex <- x$v$inftimes[hostID]
#   
#   # remove sample edges from hostID, keep parent nodes and endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   sampleparentnodes <- c()
#   for(samplenode in edgessample) {
#     parentnode <- x$v$nodeparents[samplenode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), samplenode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#     }
#     sampleparentnodes <- c(parentnode, sampleparentnodes)
#   }
#   
#   # remove sample edges from new index, keep parent nodes and endtimes
#   edgessamplenewindex <- which(x$v$nodehosts == newindex & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimesnewindex <- x$v$nodetimes[edgessamplenewindex]
#   sampleparentnodesnewindex <- c()
#   for(samplenode in edgessamplenewindex) {
#     parentnode <- x$v$nodeparents[samplenode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), samplenode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#     }
#     sampleparentnodesnewindex <- c(parentnode, sampleparentnodesnewindex)
#   }
#   
#   # change transmission tree and swap minitrees
#   infecteeshostID <- which(x$v$infectors == hostID)
#   infecteesnewindex <- which(x$v$infectors == newindex)
#   x$v$infectors[infecteeshostID] <- newindex
#   x$v$infectors[infecteesnewindex] <- hostID
#   x$v$infectors[hostID] <- newindex
#   x$v$infectors[newindex] <- 0
#   x$v$inftimes[hostID] <- newinftime_hostID
#   x$v$inftimes[newindex] <- newinftime_newindex
#   coalnodeshostID <- which(x$v$nodehosts == hostID & x$v$nodetypes == "c")
#   coalnodesnewindex <- which(x$v$nodehosts == newindex & x$v$nodetypes == "c")
#   x$v$nodehosts[coalnodeshostID] <- newindex
#   x$v$nodehosts[coalnodesnewindex] <- hostID
#   x$v$nodehosts[sampleparentnodesnewindex] <- newindex
#   x$v$nodehosts[sampleparentnodes] <- hostID
#   
#   
#   ### add sampling edges to hostID
#   # edges entering hostID, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs, 
#                                       function(y) hostID %in% c(.ptr(x$v$infectors, y))))
#   edgesend <- infecteeedges(x, hostID, transmissioncluster)
#   edgesendhosts <- x$v$nodehosts[edgesend]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgesendhosts %in% infecteeIDs))) {
#     edgesendhosts[tochange] <- x$v$infectors[edgesendhosts][tochange]
#   }
#   edgesendtimes <- x$v$inftimes[edgesendhosts]
#   
#   # coalescentnodes in currenthostID, with endtimes
#   coalescentnodes <- setdiff(which(x$v$nodehosts == hostID & x$v$nodetypes == "c"), sampleparentnodes)
#   coalescenttimes <- x$v$nodetimes[coalescentnodes]
#   
#   # new coalescenttimes
#   coalescenttimesnew <- c()
#   sampletimesadded <- c()
#   for(le in 1:length(edgessample)) {
#     coalescenttimesnew <- c(coalescenttimesnew, 
#                             x$v$inftimes[hostID] + 
#                               sample_singlecoaltime(c(edgesendtimes, sampletimesadded) - x$v$inftimes[hostID],
#                                                     c(coalescenttimes, coalescenttimesnew) - x$v$inftimes[hostID],
#                                                     edgesampletimes[le] - x$v$inftimes[hostID], x$p))
#     sampletimesadded <- c(sampletimesadded, edgesampletimes[le])
#   }
#   
#   # old within-host minitree
#   childnodes <- c(edgesend, coalescentnodes)
#   parentnodes <- x$v$nodeparents[childnodes]
#   childnodestimes <- c(edgesendtimes, coalescenttimes)
#   
#   # place new edges into minitree
#   for(le in 1:length(edgessample)) {
#     if(coalescenttimesnew[le] > x$v$inftimes[hostID]) {
#       newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
#       childnodes <- c(childnodes, edgessample[le], sampleparentnodes[le])
#       parentnodes <- c(parentnodes, sampleparentnodes[le], parentnodes[childnodes == newchildnode])
#       childnodestimes <- c(childnodestimes, edgesampletimes[le], coalescenttimesnew[le])
#       parentnodes[childnodes == newchildnode] <- sampleparentnodes[le]
#     }
#   }
#   nodestoinfector <- sampleparentnodes[coalescenttimesnew < x$v$inftimes[hostID]]
#   
#   # change phybreak object
#   x$v$nodetimes[sampleparentnodes] <- coalescenttimesnew
#   x$v$nodeparents[childnodes] <- parentnodes
#   x$v$nodehosts[nodestoinfector] <- newindex
#   x$v$nodetimes[nodestoinfector] <- x$v$inftimes[hostID]
#   x$v$nodeparents[nodestoinfector] <- -1
#   
#   
#   ### add new transmission edges and sampling edges to new index
#   # new edges entering new index
#   edgesstartnew <- which(x$v$nodeparents == -1)
#   edgesendnew <- which(x$v$nodeparents %in% edgesstartnew)
#   edgesstartnew <- x$v$nodeparents[edgesendnew]
#   
#   # old edges entering new index, with endtimes
#   transmissioncluster <- 1:x$p$obs
#   edgesend <- infecteeedges(x, newindex, transmissioncluster)
#   edgesendhosts <- x$v$nodehosts[edgesend]
#   infecteeIDs <- which(x$v$infectors == newindex)
#   while(any(tochange <- !(edgesendhosts %in% infecteeIDs))) {
#     edgesendhosts[tochange] <- x$v$infectors[edgesendhosts][tochange]
#   }
#   edgesendtimes <- x$v$inftimes[edgesendhosts]
#   
#   # coalescentnodes in new index, with endtimes
#   coalescentnodes <- setdiff(which(x$v$nodehosts == newindex & x$v$nodetypes == "c"), c(edgesstartnew, sampleparentnodesnewindex))
#   coalescenttimes <- x$v$nodetimes[coalescentnodes]
#   
#   # new coalescenttimes
#   coalescenttimesnew_sample <- c()
#   coalescenttimesnew_transmission <- c()
#   sampletimesadded <- c()
#   transmissiontimesadded <- c()
#   for(le in 1:length(edgessamplenewindex)) {
#     coalescenttimesnew_sample <- c(coalescenttimesnew_sample, 
#                                    x$v$inftimes[newindex] + 
#                                      sample_singlecoaltime(c(edgesendtimes, sampletimesadded) - x$v$inftimes[newindex],
#                                                            c(coalescenttimes, coalescenttimesnew_sample) - x$v$inftimes[newindex],
#                                                            edgesampletimesnewindex[le] - x$v$inftimes[newindex], x$p))
#     sampletimesadded <- c(sampletimesadded, edgesampletimesnewindex[le])
#   }
#   for(le in 1:length(edgesendnew)) {
#     coalescenttimesnew_transmission <- c(coalescenttimesnew_transmission, 
#                                          x$v$inftimes[newindex] + 
#                                            sample_singlecoaltime(c(edgesendtimes, sampletimesadded, transmissiontimesadded) - x$v$inftimes[newindex],
#                                                                  c(coalescenttimes, coalescenttimesnew_sample, coalescenttimesnew_transmission) - x$v$inftimes[newindex],
#                                                                  x$v$inftimes[hostID] - x$v$inftimes[newindex], x$p))
#     transmissiontimesadded <- c(transmissiontimesadded, x$v$inftimes[hostID])
#   }
#   
#   # old within-host minitree
#   childnodes <- c(edgesend, coalescentnodes)
#   parentnodes <- x$v$nodeparents[childnodes]
#   childnodestimes <- c(edgesendtimes, coalescenttimes)
#   
#   # place new edges into minitree
#   for(le in 1:length(edgessamplenewindex)) {
#     if(coalescenttimesnew_sample[le] > x$v$inftimes[newindex]) {
#       newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew_sample[le])
#       childnodes <- c(childnodes, edgessamplenewindex[le], sampleparentnodesnewindex[le])
#       parentnodes <- c(parentnodes, sampleparentnodesnewindex[le], parentnodes[childnodes == newchildnode])
#       childnodestimes <- c(childnodestimes, sampletimesadded[le], coalescenttimesnew_sample[le])
#       parentnodes[childnodes == newchildnode] <- sampleparentnodesnewindex[le]
#     }
#   }
#   for(le in 1:length(edgesendnew)) {
#     if(coalescenttimesnew_transmission[le] > x$v$inftimes[newindex]) {
#       newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew_transmission[le])
#       childnodes <- c(childnodes, edgesendnew[le], edgesstartnew[le])
#       parentnodes <- c(parentnodes, edgesstartnew[le], parentnodes[childnodes == newchildnode])
#       childnodestimes <- c(childnodestimes, transmissiontimesadded[le], coalescenttimesnew_transmission[le])
#       parentnodes[childnodes == newchildnode] <- edgesstartnew[le]
#     }
#   }
#   nodestoinfector <- sampleparentnodes[coalescenttimesnew < x$v$inftimes[newindex]]
#   
#   # change phybreak object
#   x$v$nodetimes[sampleparentnodesnewindex] <- coalescenttimesnew_sample
#   x$v$nodetimes[edgesstartnew] <- coalescenttimesnew_transmission
#   x$v$nodeparents[childnodes] <- parentnodes
#   x$v$nodehosts[nodestoinfector] <- 0
#   x$v$nodetimes[nodestoinfector] <- x$v$inftimes[newindex]
#   x$v$nodeparents[nodestoinfector] <- -1
#   
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   rewireEXP_connectindex(x, newindex)
# }  
#   
# 
# 
# ###path F2
# rewireEXP_swapcomplete <- function(x, hostID) {
#   # x is phybreakenvironment
#   
#   # identify new infector and new infection times
#   infectorID <- x$v$infector[hostID]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   newinfector <- infecteeIDs[which(x$v$inftimes[infecteeIDs] == min(x$v$inftimes[infecteeIDs]))]
#   newinftime_hostID <- x$v$inftimes[newinfector]
#   newinftime_newinfector <- x$v$inftimes[hostID]
#   
#   # remove sample edges from hostID, keep parent nodes and endtimes
#   edgessample <- which(x$v$nodehosts == hostID & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimes <- x$v$nodetimes[edgessample]
#   sampleparentnodes <- c()
#   for(samplenode in edgessample) {
#     parentnode <- x$v$nodeparents[samplenode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), samplenode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#     }
#     sampleparentnodes <- c(parentnode, sampleparentnodes)
#   }
#   
#   # remove sample edges from new infector, keep parent nodes and endtimes
#   edgessamplenewinfector <- which(x$v$nodehosts == newinfector & x$v$nodetypes %in% c("s", "x"))
#   edgesampletimesnewinfector <- x$v$nodetimes[edgessamplenewinfector]
#   sampleparentnodesnewinfector <- c()
#   for(samplenode in edgessamplenewinfector) {
#     parentnode <- x$v$nodeparents[samplenode]
#     if(parentnode != 0) {
#       second_childnode <- setdiff(which(x$v$nodeparents == parentnode), samplenode)
#       x$v$nodeparents[second_childnode] <- x$v$nodeparents[parentnode]
#     }
#     sampleparentnodesnewinfector <- c(parentnode, sampleparentnodesnewinfector)
#   }
#   
#   # change transmission tree and swap minitrees
#   infecteeshostID <- which(x$v$infectors == hostID)
#   infecteesnewinfector <- which(x$v$infectors == newinfector)
#   x$v$infectors[infecteeshostID] <- newinfector
#   x$v$infectors[infecteesnewinfector] <- hostID
#   x$v$infectors[hostID] <- newinfector
#   x$v$infectors[newinfector] <- infectorID
#   x$v$inftimes[hostID] <- newinftime_hostID
#   x$v$inftimes[newinfector] <- newinftime_newinfector
#   coalnodeshostID <- which(x$v$nodehosts == hostID & x$v$nodetypes == "c")
#   coalnodesnewinfector <- which(x$v$nodehosts == newinfector & x$v$nodetypes == "c")
#   x$v$nodehosts[coalnodeshostID] <- newinfector
#   x$v$nodehosts[coalnodesnewinfector] <- hostID
#   x$v$nodehosts[sampleparentnodesnewinfector] <- newinfector
#   x$v$nodehosts[sampleparentnodes] <- hostID
#   
#   
#   ### add sampling edges to hostID
#   # edges entering hostID, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs, 
#                                       function(y) hostID %in% c(.ptr(x$v$infectors, y))))
#   edgesend <- infecteeedges(x, hostID, transmissioncluster)
#   edgesendhosts <- x$v$nodehosts[edgesend]
#   infecteeIDs <- which(x$v$infectors == hostID)
#   while(any(tochange <- !(edgesendhosts %in% infecteeIDs))) {
#     edgesendhosts[tochange] <- x$v$infectors[edgesendhosts][tochange]
#   }
#   edgesendtimes <- x$v$inftimes[edgesendhosts]
#   
#   # coalescentnodes in currenthostID, with endtimes
#   coalescentnodes <- setdiff(which(x$v$nodehosts == hostID & x$v$nodetypes == "c"), sampleparentnodes)
#   coalescenttimes <- x$v$nodetimes[coalescentnodes]
#   
#   # new coalescenttimes
#   coalescenttimesnew <- c()
#   sampletimesadded <- c()
#   for(le in 1:length(edgessample)) {
#     coalescenttimesnew <- c(coalescenttimesnew, 
#                             x$v$inftimes[hostID] + 
#                               sample_singlecoaltime(c(edgesendtimes, sampletimesadded) - x$v$inftimes[hostID],
#                                                     c(coalescenttimes, coalescenttimesnew) - x$v$inftimes[hostID],
#                                                     edgesampletimes[le] - x$v$inftimes[hostID], x$p))
#     sampletimesadded <- c(sampletimesadded, edgesampletimes[le])
#   }
#   
#   # old within-host minitree
#   childnodes <- c(edgesend, coalescentnodes)
#   parentnodes <- x$v$nodeparents[childnodes]
#   childnodestimes <- c(edgesendtimes, coalescenttimes)
#   
#   # place new edges into minitree
#   for(le in 1:length(edgessample)) {
#     if(coalescenttimesnew[le] > x$v$inftimes[hostID]) {
#       newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew[le])
#       childnodes <- c(childnodes, edgessample[le], sampleparentnodes[le])
#       parentnodes <- c(parentnodes, sampleparentnodes[le], parentnodes[childnodes == newchildnode])
#       childnodestimes <- c(childnodestimes, edgesampletimes[le], coalescenttimesnew[le])
#       parentnodes[childnodes == newchildnode] <- sampleparentnodes[le]
#     }
#   }
#   nodestoinfector <- sampleparentnodes[coalescenttimesnew < x$v$inftimes[hostID]]
#   
#   # change phybreak object
#   x$v$nodetimes[sampleparentnodes] <- coalescenttimesnew
#   x$v$nodeparents[childnodes] <- parentnodes
#   x$v$nodehosts[nodestoinfector] <- newinfector
#   x$v$nodetimes[nodestoinfector] <- x$v$inftimes[hostID]
#   x$v$nodeparents[nodestoinfector] <- -1
#   
#   
#   ### add transmission and sampling edges to new infector
#   # new edges entering new infector
#   edgesstartnew <- which(x$v$nodeparents == -1)
#   edgesendnew <- which(x$v$nodeparents %in% edgesstartnew)
#   edgesstartnew <- x$v$nodeparents[edgesendnew]
#   
#   # old edges entering new index, with endtimes
#   transmissioncluster <- which(sapply(1:x$p$obs, 
#                                       function(y) newinfector %in% c(.ptr(x$v$infectors, y))))
#   edgesend <- setdiff(infecteeedges(x, newinfector, transmissioncluster), edgesendnew)
#   edgesendhosts <- x$v$nodehosts[edgesend]
#   infecteeIDs <- which(x$v$infectors == newinfector)
#   while(any(tochange <- !(edgesendhosts %in% infecteeIDs))) {
#     edgesendhosts[tochange] <- x$v$infectors[edgesendhosts][tochange]
#   }
#   edgesendtimes <- x$v$inftimes[edgesendhosts]
#   
#   # coalescentnodes in currenthostID, with endtimes
#   coalescentnodes <- setdiff(which(x$v$nodehosts == newinfector & x$v$nodetypes == "c"), c(edgesstartnew, sampleparentnodesnewinfector))
#   coalescenttimes <- x$v$nodetimes[coalescentnodes]
#   
#   # new coalescenttimes
#   coalescenttimesnew_sample <- c()
#   coalescenttimesnew_transmission <- c()
#   sampletimesadded <- c()
#   transmissiontimesadded <- c()
#   for(le in 1:length(edgessamplenewinfector)) {
#     coalescenttimesnew_sample <- c(coalescenttimesnew_sample, 
#                                    x$v$inftimes[newinfector] + 
#                                      sample_singlecoaltime(c(edgesendtimes, sampletimesadded) - x$v$inftimes[newinfector],
#                                                            c(coalescenttimes, coalescenttimesnew_sample) - x$v$inftimes[newinfector],
#                                                            edgesampletimesnewinfector[le] - x$v$inftimes[newinfector], x$p))
#     sampletimesadded <- c(sampletimesadded, edgesampletimesnewinfector[le])
#   }
#   if(length(edgesendnew) > 0) {
#     for(le in 1:length(edgesendnew)) {
#       coalescenttimesnew_transmission <- c(coalescenttimesnew_transmission, 
#                                            x$v$inftimes[newinfector] + 
#                                              sample_singlecoaltime(c(edgesendtimes, sampletimesadded, transmissiontimesadded) - x$v$inftimes[newinfector],
#                                                                    c(coalescenttimes, coalescenttimesnew_sample, coalescenttimesnew_transmission) - x$v$inftimes[newinfector],
#                                                                    x$v$inftimes[hostID] - x$v$inftimes[newinfector], x$p))
#       transmissiontimesadded <- c(transmissiontimesadded, x$v$inftimes[hostID])
#     }
#   }
#   
#   # old within-host minitree
#   childnodes <- c(edgesend, coalescentnodes)
#   parentnodes <- x$v$nodeparents[childnodes]
#   childnodestimes <- c(edgesendtimes, coalescenttimes)
#   
#   # place new edges into minitree
#   for(le in 1:length(edgessamplenewinfector)) {
#     if(coalescenttimesnew_sample[le] > x$v$inftimes[newinfector]) {
#       newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew_sample[le])
#       childnodes <- c(childnodes, edgessamplenewinfector[le], sampleparentnodesnewinfector[le])
#       parentnodes <- c(parentnodes, sampleparentnodesnewinfector[le], parentnodes[childnodes == newchildnode])
#       childnodestimes <- c(childnodestimes, sampletimesadded[le], coalescenttimesnew_sample[le])
#       parentnodes[childnodes == newchildnode] <- sampleparentnodesnewinfector[le]
#     }
#   }
#   if(length(edgesendnew) > 0) {
#     for(le in 1:length(edgesendnew)) {
#       if(coalescenttimesnew_transmission[le] > x$v$inftimes[newinfector]) {
#         newchildnode <- sample_singlechildnode(childnodes, parentnodes, childnodestimes, coalescenttimesnew_transmission[le])
#         childnodes <- c(childnodes, edgesendnew[le], edgesstartnew[le])
#         parentnodes <- c(parentnodes, edgesstartnew[le], parentnodes[childnodes == newchildnode])
#         childnodestimes <- c(childnodestimes, transmissiontimesadded[le], coalescenttimesnew_transmission[le])
#         parentnodes[childnodes == newchildnode] <- edgesstartnew[le]
#       }
#     }
#   }
#   nodestoinfector <- c(sampleparentnodes[coalescenttimesnew_sample < x$v$inftimes[newinfector]],
#                        sampleparentnodes[coalescenttimesnew_transmission < x$v$inftimes[newinfector]])
#   
#   # change phybreak object
#   x$v$nodetimes[sampleparentnodesnewinfector] <- coalescenttimesnew_sample
#   x$v$nodetimes[edgesstartnew] <- coalescenttimesnew_transmission
#   x$v$nodeparents[childnodes] <- parentnodes
#   x$v$nodehosts[nodestoinfector] <- infectorID
#   x$v$nodetimes[nodestoinfector] <- x$v$inftimes[newinfector]
#   x$v$nodeparents[nodestoinfector] <- -1
#   
#   # connect edges in the bottleneck to the phylogenetic tree in the infector
#   rewireEXP_connectinfector(x, newinfector)
# }  
# 
# 
# 
# 
