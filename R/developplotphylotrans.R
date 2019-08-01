# library(phybreak)
# simul <- sim_phybreak(10, samplesperhost = 5, wh.level = 1, wh.bottleneck = "wide")
# 
# plotje <- plotPhyloTrans(simul)
# 
# plotTrans(simul)
# 
# hostorder <- c(10,8,7,9,5,6,4,2,3,1,0)
# 
# treelist <- list()
# for(i in 1:length(plotje$v$nodeparents)) {
#   if(plotje$v$nodetypes[i] != "c") {
#     treelist <- c(
#       treelist,
#       list(
#         list(
#           x1vec = plotje$v$nodetimes[plotje$v$nodeparents[i]],
#           x2vec = plotje$v$nodetimes[i],
#           y1vec = 0,
#           y2vec = 0,
#           nodevec = i,
#           rootnode = plotje$v$nodeparents[i],
#           xroot = plotje$v$nodetimes[plotje$v$nodeparents[i]],
#           yscore = if(plotje$v$nodetypes[i] %in% c("t", "b")) {
#             if(which(hostorder == plotje$v$nodehosts[i]) >
#                which(hostorder == plotje$v$nodehosts[
#                  which(plotje$v$nodeparents == i)
#                  ])) {
#               -1
#             } else 1
#           } else 0,
#           host = plotje$v$nodehosts[i],
#           completeYN = if(plotje$v$nodetypes[plotje$v$nodeparents[i]] %in% c("t", "b")) {
#             TRUE
#           } else FALSE
#         )
#       )
#     )
#   }
# }
# hosttrees <- list()
# while(length(treelist) > 0) {
#   xroots <- sapply(treelist, function(x) x$xroot)
#   tojoin <- which(xroots == max(xroots))
#   
#   trees2join <- treelist[tojoin]
#   treelist <- treelist[-tojoin]
#   
#   #not yet ordered by *time* of secondary transmission
#   treeorder <- order(sapply(trees2join, function(x) x$yscore))
#   trees2join <- trees2join[treeorder]
#   
#   for(i in 1:length(trees2join)) {
#     trees2join$y1vec <- trees2join$y1vec
#   }
#   
#   newtree <- list(
#     x1vec = c(sapply(trees2join, function(x) x$x1vec)),
#     x2vec = c(sapply(trees2join, function(x) x$x2vec)),
#     y1vec = c(sapply(trees2join, function(x) x$y1vec)),
#     y2vec = c(sapply(trees2join, function(x) x$y2vec)),
#     nodevec = c(sapply(trees2join, function(x) x$nodevec))
#   )
# }
# xroots <- sapply(treelist, function(x) x$xroot)
