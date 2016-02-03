# phybreak.env <- new.env()
#
#
#
# ###For likelihood calculation
#
# #helper functions
# .makelik.all <- function(SNPs,obs) {
#   SNPnrs <- SNPs == 0
#   SNPnrs[SNPs=="a"] <- 1
#   SNPnrs[SNPs=="c"] <- 2
#   SNPnrs[SNPs=="g"] <- 3
#   SNPnrs[SNPs=="t"] <- 4
#   phybreak.env$.lik.all <- array(NA, dim = c(3*obs-1,ncol(SNPs),4))
#   phybreak.env$.lik.all[1:obs,,] <- 0
#   phybreak.env$.lik.all[cbind(rep(1:obs,ncol(SNPs)),
#                               rep(1:ncol(SNPs),each=obs),
#                               as.vector(SNPnrs))] <- 1
# }
# .clearlik.all <- function(hostID, obs, nodehosts) {
#   if(hostID > 0 & !is.na(phybreak.env$.lik.all[hostID + 2*obs -1, 1, 1])) {
#     phybreak.env$.lik.all[hostID + 2*obs - 1,1,1] <- NA
#     .clearlik.all(nodehosts[hostID + 2*obs - 1], obs, nodehosts)
#   }
# }
# .liknucl <- function(off, le, mutrate) {
#   (1/4)*rowSums(off)+(off-(1/4)*rowSums(off))*exp(-mutrate*le)
# }
# .likseqnode <- function(nID, mutrate, wr) {
#   if(phybreak.env$.n.ty[nID] == "c") {
#     iees <- which(phybreak.env$.n.pa == nID)
#     newlik <- .liknucl(
#       .likseqnode(iees[1], mutrate, wr),
#       phybreak.env$.n.ti[iees[1]] - phybreak.env$.n.ti[nID],
#       mutrate
#     ) *
#       .liknucl(
#         .likseqnode(iees[2], mutrate, wr),
#         phybreak.env$.n.ti[iees[2]] - phybreak.env$.n.ti[nID],
#         mutrate
#       )
#     return(newlik)
#   }
#   if(is.na(phybreak.env$.lik.all[nID, 1, 1])) {
#     iees <- which(phybreak.env$.n.pa == nID)
#     newlik <- .liknucl(
#       .likseqnode(iees, mutrate, wr),
#       phybreak.env$.n.ti[iees] - phybreak.env$.n.ti[nID],
#       mutrate
#     )
#     if(wr) phybreak.env$.lik.all[nID,,] <- newlik
#     return(newlik)
#   }
#   return(phybreak.env$.lik.all[nID,,])
# }
#
# #calculating likelihood with current parameters
# .lik.sequences <- function(mutrate, SNPs, SNPfreqs, obs, nodetimes, nodehosts, nodeparents, nodetypes,
#                            withinroutine = FALSE) {
#   phybreak.env$.n.ti <- nodetimes
#   phybreak.env$.n.pa <- nodeparents
#   phybreak.env$.n.ty <- nodetypes
#   if(!withinroutine | !exists(".lik.all", envir = phybreak.env)) {
#     if(exists(".lik.all", envir = phybreak.env)) likallcp <- phybreak.env$.lik.all
#     .makelik.all(SNPs,obs)
#     ans <- sum(log(.25*rowSums(.likseqnode(which(nodehosts==0), mutrate, wr=FALSE)) ) * SNPfreqs)
#     if(exists("likallcp")) phybreak.env$.lik.all <- likallcp else rm(.lik.all,envir = phybreak.env)
#   }
#   else ans <- sum(log(.25*rowSums(.likseqnode(which(nodehosts==0), mutrate, wr=TRUE)) ) * SNPfreqs)
#   rm(.n.ti,.n.pa,.n.ty, envir=phybreak.env)
#   return(ans)
# }
#
#
# #the actual likelihood
# lik.phybreak <- function(phybreak.object, seq = TRUE, t.g = FALSE, t.s = FALSE, withinroutine = FALSE) {
#   (if(seq) with(phybreak.object,.lik.sequences(p$mu,d$SNP,d$SNPfr,p$obs,v$nodetimes,
#                                                v$nodehosts,v$nodeparents,v$nodetypes,withinroutine)) else 0) +
#     (if(t.g) with(phybreak.object,.lik.gentimes(p$obs,p$shape.gen,p$mean.gen,v$nodetimes,v$nodehosts,v$nodetypes)) else 0) +
#     (if(t.s) with(phybreak.object,.lik.sampletimes(p$shape.sample,p$mean.sample,v$nodetimes,v$nodetypes)) else 0)
# }
#
#
#
#
# ###For updating the global outbreak variables
# burnin.phybreak <- function(phybreak.object, nburnin) {
#   parked.samples <- phybreak.object$s
#   .makelik.all(phybreak.object$d$SNP, phybreak.object$p$obs)
#   res <- phybreak.object
#   res$s <- c()
#
#   tdifpropose <<- 0
#   tdiflikcalc <<- 0
#   tdifhosts <- 0
#   tdifrest <- 0
#
#   for(rep in 1:nburnin) {
#     begi <- Sys.time()
#     for(i in sample(phybreak.object$p$obs)) {
#       res <- .updatehost(i, res)
#     }
#     tdifhosts <- tdifhosts + Sys.time() - begi
#     begi <- Sys.time()
#     if(res$h$est.mG) res <- .update.mG(res)
#     if(res$h$est.mS) res <- .update.mS(res)
#     res <- .update.mu(res)
#     tdifrest <- tdifrest + Sys.time() - begi
#   }
#   rm(.lik.all, envir = phybreak.env)
#
#   res$s <- parked.samples
#
#   units(tdifhosts) <- "secs"
#   units(tdifrest) <- "secs"
#   units(tdifpropose) <- "secs"
#   units(tdiflikcalc) <- "secs"
#
#   print(c("host updates: ", tdifhosts))
#   print(c("rest updates: ", tdifrest))
#   print(c("proposing: ", tdifpropose))
#   print(c("likelihood calculation: ", tdiflikcalc))
#
#   return(res)
# }
#
# sample.phybreak <- function(phybreak.object, nsample, thin) {
#   s.post <- list(
#     nodetimes = with(phybreak.object,cbind(s$nodetimes, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
#     nodehosts = with(phybreak.object,cbind(s$nodehosts, matrix(NA, nrow=2*p$obs - 1, ncol=nsample))),
#     nodeparents = with(phybreak.object,cbind(s$nodeparents, matrix(NA, nrow=3*p$obs - 1, ncol=nsample))),
#     mu = c(phybreak.object$s$mu, rep(NA, nsample)),
#     mG = c(phybreak.object$s$mG, rep(NA, nsample)),
#     mS = c(phybreak.object$s$mS, rep(NA, nsample)),
#     logLik = c(phybreak.object$s$logLik, rep(NA, nsample))
#   )
#   res <- phybreak.object
#   res$s <- c()
#   .makelik.all(phybreak.object$d$SNP, phybreak.object$p$obs)
#
#   for(sa in tail(1:length(s.post$mu), nsample)) {
#     for(rep in 1:thin) {
#       for(i in sample(phybreak.object$p$obs)) {
#         res <- .updatehost(i, res)
#       }
#       if(res$h$est.mG) res <- .update.mG(res)
#       if(res$h$est.mS) res <- .update.mS(res)
#       res <- .update.mu(res)
#     }
#     s.post$nodetimes[, sa] <- tail(res$v$nodetimes, -res$p$obs)
#     s.post$nodehosts[, sa] <- tail(res$v$nodehosts, -res$p$obs)
#     s.post$nodeparents[, sa] <- res$v$nodeparents
#     s.post$mu[sa] <- res$p$mu
#     s.post$mG[sa] <- res$p$mean.gen
#     s.post$mS[sa] <- res$p$mean.sample
#     s.post$logLik[sa] <- lik.phybreak(res, seq = TRUE, t.g = TRUE, t.s = TRUE, withinroutine = TRUE)
#   }
#
#   rm(.lik.all, envir = phybreak.env)
#
#   res$s <- s.post
#
#   return(res)
#
# }
#
#
#
# .updatehost <- function(hostID, phybreak.object) {
#   with(phybreak.object, {
#     #    tinf.prop <- v$nodetimes[hostID] - rgamma(1,shape=p$shape.sample,scale=p$mean.sample/p$shape.sample)
#     tinf.prop <- v$nodetimes[hostID] - rgamma(1,shape=2,scale=p$mean.sample/2)
#
#     ##going down the decision tree
#     if(tinf.prop < suppressWarnings(min(v$nodetimes[v$nodetypes == "t" & v$nodehosts == hostID]))) {
#       #Q1=Y : tinf.prop before first infectee
#       if(v$nodehosts[hostID + 2*p$obs - 1] == 0) {
#         #Q12=YY : tinf.prop before first infectee & hostID is index
#         return(.updatepathA(hostID, tinf.prop, phybreak.object))
#       } else {
#         #Q12=YN : tinf.prop before first infectee & hostID is not index
#         if(tinf.prop < v$nodetimes[v$nodehosts == 0]) {
#           #Q123=YNY : hostID is not index & tinf.prop is before infection of index
#           iees <- which(v$nodehosts == (which(v$nodehosts == 0) - 2*p$obs + 1) & v$nodetypes == "t")
#           if(hostID == iees[order(v$nodetimes[iees])][1] - 2*p$obs + 1) {
#             #Q1234=YNYY : tinf.prop before infection of index & hostID is index's first infectee
#             return(.updatepathB(hostID, tinf.prop, phybreak.object))
#           } else {
#             #Q1234=YNYN : tinf.prop before infection of index & hostID is not index nor index's first infectee
#             ### DO NOTHING
#             return(phybreak.object)
#           }
#         } else {
#           #Q123=YNN : tinf.prop before first infectee and after infection of index
#           return(.updatepathC(hostID, tinf.prop, phybreak.object))
#         }
#       }
#     } else
#     {
#       #Q1=N : tinf.prop after first infectee
#       if(v$nodehosts[hostID + 2*p$obs - 1] == 0) {
#         #Q12=NY : tinf.prop after first infectee & hostID is index
#         if(tinf.prop > sort(tail(v$nodetimes,p$obs))[3]) {
#           #Q123=NYY : hostID is index & tinf.prop is after third infection in outbreak
#           ### DO NOTHING
#           return(phybreak.object)
#         } else {
#           #Q123=NYN : hostID is index & tinf.prop is after first infectee but before third infection in outbreak
#           return(.updatepathD(hostID, tinf.prop, phybreak.object))
#         }
#       } else {
#         #Q12=NN : tinf.prop after first infectee & hostID is not index
#         return(.updatepathE(hostID, phybreak.object))
#       }
#     }
#   })
# }
#
# #helper functions for updating variables
# {
# .updatepathA <- function(hostID, tinf.prop, phybreak.object) {
#   ###save current state
#   lik.all.current <- phybreak.env$.lik.all
#   lik.current <- lik.phybreak(phybreak.object, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE)
#   ###change to proposal state
#
#   beg <- Sys.time()
#
#   phybreak.obj.prop <- within(phybreak.object,{
#     v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
#
#     v$nodetimes[v$nodehosts == hostID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
#       v$nodetimes[hostID + 2*p$obs - 1] +                      #...to the infection time +
#       .samplecoaltimes(v$nodetimes[v$nodehosts == hostID & v$nodetypes != "c"] - v$nodetimes[hostID + 2*p$obs - 1],
#                        p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times
#
#     v$nodeparents[v$nodehosts == hostID] <-     #change the parent nodes of all nodes in hostID...
#       .sampletopology(which(v$nodehosts == hostID), v$nodetimes[v$nodehosts == hostID], v$nodetypes[v$nodehosts == hostID],
#                       hostID + 2*p$obs - 1, p$wh.model)
#     #...to a correct topology, randomized where possible
#
#     .clearlik.all(hostID, p$obs, v$nodehosts)
#   })
#
#   tdifpropose <<- tdifpropose + Sys.time() - beg
#   beg <- Sys.time()
#
#
#   ###calculate acceptance probability
#   logproposalratio <-     with(phybreak.object,
#                                dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                       shape = 2, scale = p$mean.sample/2,log=TRUE)) -
#     with(phybreak.obj.prop, dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                    shape = 2, scale = p$mean.sample/2, log=TRUE))
#
#   logaccprob <- lik.phybreak(phybreak.obj.prop, seq=TRUE,
#                              t.g=TRUE, t.s=TRUE, withinroutine = TRUE) +
#     logproposalratio - lik.current
#
#   tdiflikcalc <<- tdiflikcalc + Sys.time() - beg
#
#
#   if(runif(1,0,1) > exp(logaccprob)) {
#     phybreak.env$.lik.all <- lik.all.current
#     return(phybreak.object)
#   } else return(phybreak.obj.prop)
#
# }
#
# .updatepathB <- function(hostID, tinf.prop, phybreak.object) {
#   ###save current state
#   lik.all.current <- phybreak.env$.lik.all
#   lik.current <- lik.phybreak(phybreak.object, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE)
#
#
#   beg <- Sys.time()
#
#
#   ###change to proposal state
#   index.current.ID <- with(phybreak.object, {
#     which(v$nodehosts == 0) - 2*p$obs + 1
#   })
#
#
#   phybreak.obj.prop <- within(phybreak.object,{
#     v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
#
#     v$nodehosts[hostID + 2*p$obs - 1] <- 0
#     v$nodeparents[hostID + 2*p$obs - 1] <- 0
#     v$nodehosts[index.current.ID + 2*p$obs - 1] <- hostID
#     v$nodehosts[v$nodehosts == index.current.ID & v$nodetypes == "c"][1] <- hostID
#
#     for(ID in c(hostID,index.current.ID)) {
#       v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
#         v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
#         .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
#                          p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times
#
#       v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
#         .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
#                         ID + 2*p$obs - 1, p$wh.model)
#       #...to a correct topology, randomized where possible
#
#     }
#     rm(ID)
#
#   })
#
#   tdifpropose <<- tdifpropose + Sys.time() - beg
#   beg <- Sys.time()
#
#   ###calculate acceptance probability
#   with(phybreak.object, {
#     .clearlik.all(hostID, p$obs, v$nodehosts)
#     .clearlik.all(index.current.ID, p$obs, v$nodehosts)
#
#   })
#   logproposalratio <- with(phybreak.object,
#                            dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                   shape = 2, scale = p$mean.sample/2,log=TRUE)) -
#     with(phybreak.obj.prop, dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                    shape = 2, scale = p$mean.sample/2, log=TRUE))
#
#
#
#   logaccprob <- lik.phybreak(phybreak.obj.prop, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE) +
#     logproposalratio - lik.current
#
#   tdiflikcalc <<- tdiflikcalc + Sys.time() - beg
#
#   if(runif(1,0,1) > exp(logaccprob)) {
#     phybreak.env$.lik.all <- lik.all.current
#     return(phybreak.object)
#   } else return(phybreak.obj.prop)
#
# }
#
# .updatepathC <- function(hostID, tinf.prop, phybreak.object) {
#   ###save current state
#   lik.all.current <- phybreak.env$.lik.all
#   lik.current <- lik.phybreak(phybreak.object, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE)
#
#   beg <- Sys.time()
#
#
#   ###change to proposal state
#   infector.current.ID <- with(phybreak.object, v$nodehosts[hostID + 2*p$obs - 1])
#   dens.infectorproposal <- with(phybreak.object,
#                                 (dgamma(tinf.prop - tail(v$nodetimes,p$obs),
#                                         shape = p$shape.gen,
#                                         scale = p$mean.gen/p$shape.gen)) +
#                                   (tinf.prop - tail(v$nodetimes,p$obs) > 0) / h$dist[hostID,])
#   dens.infectorproposal[hostID] <- 0
#   infector.proposed.ID <- sample(phybreak.object$p$obs, 1, prob = dens.infectorproposal)
#
#
#   phybreak.obj.prop <- within(phybreak.object,{
#     v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
#     v$nodehosts[hostID + 2*p$obs - 1] <- infector.proposed.ID
#     v$nodehosts[v$nodehosts == infector.current.ID & v$nodetypes == "c"][1] <- infector.proposed.ID
#
#
#     for(ID in c(hostID,infector.current.ID,infector.proposed.ID)) {
#       v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
#         v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
#         .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
#                          p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times
#
#       v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
#         .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
#                         ID + 2*p$obs - 1, p$wh.model)
#       #...to a correct topology, randomized where possible
#
#     }
#     rm(ID)
#
#   })
#
#   tdifpropose <<- tdifpropose + Sys.time() - beg
#   beg <- Sys.time()
#
#
#   ###calculate acceptance probability
#   with(phybreak.object,{
#     .clearlik.all(hostID, p$obs, v$nodehosts)
#     .clearlik.all(infector.current.ID, p$obs, v$nodehosts)
#     .clearlik.all(infector.proposed.ID, p$obs, v$nodehosts)
#   })
#   dens.infectorcurrent <- with(phybreak.object,
#                                (dgamma(v$nodetimes[hostID + 2*p$obs - 1] -
#                                          tail(v$nodetimes,p$obs),
#                                        shape = p$shape.gen, scale = p$mean.gen/p$shape.gen)) +
#                                  (v$nodetimes[hostID + 2*p$obs - 1] - tail(v$nodetimes,p$obs) > 0) / h$dist[hostID,])
#   dens.infectorcurrent[hostID] <- 0
#   logproposalratio <- log(dens.infectorcurrent[infector.current.ID] * sum(dens.infectorproposal) /
#                             (dens.infectorproposal[infector.proposed.ID] * sum(dens.infectorcurrent))) +
#     with(phybreak.object,
#          dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                 shape = 2, scale = p$mean.sample/2,log=TRUE)) -
#     with(phybreak.obj.prop, dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                    shape = 2, scale = p$mean.sample/2, log=TRUE))
#
#   logaccprob <- lik.phybreak(phybreak.obj.prop, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE) -
#     lik.current + logproposalratio
#
#   tdiflikcalc <<- tdiflikcalc + Sys.time() - beg
#
#   if(runif(1,0,1) > exp(logaccprob)) {
#     phybreak.env$.lik.all <- lik.all.current
#     return(phybreak.object)
#   } else return(phybreak.obj.prop)
#
# }
#
# .updatepathD <- function(hostID, tinf.prop, phybreak.object) {
#   ###save current state
#   lik.all.current <- phybreak.env$.lik.all
#   lik.current <- lik.phybreak(phybreak.object, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE)
#
#
#   beg <- Sys.time()
#
#   ###change to proposal state
#
#   iees.nodeIDs <- with(phybreak.object, which(v$nodehosts == hostID & v$nodetypes == "t"))
#   infectee.current.ID <- with(phybreak.object,
#                               iees.nodeIDs[v$nodetimes[iees.nodeIDs] == min(v$nodetimes[iees.nodeIDs])] - 2*p$obs + 1)
#
#   phybreak.obj.prop <- within(phybreak.object, {
#     v$nodetimes[hostID + 2*p$obs - 1] <- tinf.prop
#
#     v$nodehosts[infectee.current.ID + 2*p$obs - 1] <- 0
#     v$nodeparents[infectee.current.ID + 2*p$obs - 1] <- 0
#     v$nodehosts[hostID + 2*p$obs - 1] <- infectee.current.ID
#     v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infectee.current.ID
#
#     for(ID in c(hostID,infectee.current.ID)) {
#       v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
#         v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
#         .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
#                          p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times
#
#       v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
#         .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
#                         ID + 2*p$obs - 1, p$wh.model)
#       #...to a correct topology, randomized where possible
#
#     }
#     rm(ID)
#
#   })
#
#   tdifpropose <<- tdifpropose + Sys.time() - beg
#   beg <- Sys.time()
#
#   ###calculate acceptance probability
#   with(phybreak.object,{
#     .clearlik.all(hostID, p$obs, v$nodehosts)
#     .clearlik.all(infectee.current.ID, p$obs, v$nodehosts)
#   })
#   logproposalratio <- with(phybreak.object,
#                            dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                   shape = 2, scale = p$mean.sample/2,log=TRUE)) -
#     with(phybreak.obj.prop, dgamma(v$nodetimes[hostID] - v$nodetimes[hostID + 2*p$obs - 1],
#                                    shape = 2, scale = p$mean.sample/2, log=TRUE))
#
#
#   logaccprob <- lik.phybreak(phybreak.obj.prop, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE) -
#     lik.current + logproposalratio
#
#   tdiflikcalc <<- tdiflikcalc + Sys.time() - beg
#
#   if(runif(1,0,1) > exp(logaccprob)) {
#     phybreak.env$.lik.all <- lik.all.current
#     return(phybreak.object)
#   } else return(phybreak.obj.prop)
#
# }
#
# .updatepathE <- function(hostID, phybreak.object) {
#   ###save current state
#   lik.all.current <- phybreak.env$.lik.all
#   lik.current <- lik.phybreak(phybreak.object, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE)
#
#
#   beg <- Sys.time()
#
#   ###change to proposal state
#   iees.nodeIDs <- with(phybreak.object, which(v$nodehosts == hostID & v$nodetypes == "t"))
#   infectee.current.ID <- with(phybreak.object,
#                               iees.nodeIDs[v$nodetimes[iees.nodeIDs] == min(v$nodetimes[iees.nodeIDs])] - 2*p$obs + 1)
#   infector.current.ID <- with(phybreak.object, v$nodehosts[hostID + 2*p$obs - 1])
#   inftime.current.hostID <- with(phybreak.object, v$nodetimes[hostID + 2*p$obs - 1])
#
#   phybreak.obj.prop <- within(phybreak.object, {
#     v$nodetimes[hostID + 2*p$obs - 1] <- v$nodetimes[infectee.current.ID + 2*p$obs - 1]
#     v$nodetimes[infectee.current.ID + 2*p$obs - 1] <- inftime.current.hostID
#
#     v$nodehosts[hostID + 2*p$obs - 1] <- infectee.current.ID
#     v$nodehosts[infectee.current.ID + 2*p$obs - 1] <- infector.current.ID
#     v$nodeparents[infectee.current.ID + 2*p$obs - 1] <- v$nodeparents[hostID + 2*p$obs - 1]
#     v$nodehosts[v$nodehosts == hostID & v$nodetypes == "c"][1] <- infectee.current.ID
#
#     for(ID in c(hostID,infectee.current.ID)) {
#       v$nodetimes[v$nodehosts == ID & v$nodetypes == "c"] <-   #change the times of the coalescence nodes in hostID...
#         v$nodetimes[ID + 2*p$obs - 1] +                      #...to the infection time +
#         .samplecoaltimes(v$nodetimes[v$nodehosts == ID & v$nodetypes != "c"] - v$nodetimes[ID + 2*p$obs - 1],
#                          p$wh.model, p$wh.lambda, p$wh.rate0)  #...sampled coalescence times
#
#       v$nodeparents[v$nodehosts == ID] <-     #change the parent nodes of all nodes in hostID...
#         .sampletopology(which(v$nodehosts == ID), v$nodetimes[v$nodehosts == ID], v$nodetypes[v$nodehosts == ID],
#                         ID + 2*p$obs - 1, p$wh.model)
#       #...to a correct topology, randomized where possible
#
#     }
#     rm(ID)
#
#   })
#
#   tdifpropose <<- tdifpropose + Sys.time() - beg
#   beg <- Sys.time()
#
#
#   ###calculate acceptance probability
#   with(phybreak.object,{
#     .clearlik.all(hostID, p$obs, v$nodehosts)
#     .clearlik.all(infectee.current.ID, p$obs, v$nodehosts)
#   })
#
#
#   logproposalratio <- with(phybreak.object,
#                            pgamma(v$nodetimes[infectee.current.ID] - v$nodetimes[infectee.current.ID + 2*p$obs - 1],
#                                   shape = 2, scale = p$mean.sample/2, log.p = TRUE) -
#                              pgamma(v$nodetimes[hostID] - v$nodetimes[infector.current.ID + 2*p$obs - 1],
#                                     shape = 2, scale = p$mean.sample/2, log.p = TRUE))
#
#
#   logaccprob <- lik.phybreak(phybreak.obj.prop, seq=TRUE, t.g=TRUE, t.s=TRUE, withinroutine = TRUE) -
#     lik.current + logproposalratio
#
#   tdiflikcalc <<- tdiflikcalc + Sys.time() - beg
#
#   if(runif(1,0,1) > exp(logaccprob)) {
#     phybreak.env$.lik.all <- lik.all.current
#     return(phybreak.object)
#   } else return(phybreak.obj.prop)
#
#
# }
#
# }
#
# ###For updating the parameters
# .update.mu <- function(phybreak.object) {
#   lik.all.current <- phybreak.env$.lik.all
#   lik.current <- lik.phybreak(phybreak.object, seq=TRUE, t.g=TRUE, t.s=FALSE, withinroutine = TRUE)
#
#   phybreak.obj.prop <- within(phybreak.object,p$mu<-exp(log(p$mu) +
#                                                           rnorm(1,0,sd(log(h$si),na.rm=TRUE))))
#
#   with(phybreak.object,phybreak.env$.lik.all[(2*p$obs):(3*p$obs-1),1,1] <- NA)
#   logaccprob <- lik.phybreak(phybreak.obj.prop, seq=TRUE, t.g=TRUE, t.s=FALSE, withinroutine = TRUE) - lik.current
#
#   if(runif(1,0,1) > exp(logaccprob)) {
#     phybreak.env$.lik.all <- lik.all.current
#     return(within(phybreak.object,h$si <- c(h$si[-1], p$mu)))
#   } else return(within(phybreak.obj.prop,h$si <- c(h$si[-1], p$mu)))
# }
#
# .update.mS <- function(phybreak.object, prior.sh = 0, prior.sc = 0) {
#   sumst <- sum(with(phybreak.object, v$nodetimes[v$nodetypes == "s"] -
#                       v$nodetimes[v$nodetypes == "t"]))
#   newmeanS <- with(phybreak.object,
#                    p$shape.sample / rgamma(1,shape = p$shape.sample * p$obs + prior.sh,
#                                            rate = sumst + prior.sc/p$shape.sample))
#   return(within(phybreak.object,
#                 p$mean.sample <- newmeanS))
#
# }
#
# .update.mG <- function(phybreak.object, prior.sh = 0, prior.sc = 0) {
#   sumgt <- sum(with(phybreak.object,
#                     v$nodetimes[v$nodetypes == "t" & v$nodehosts != 0] -
#                       v$nodetimes[v$nodetypes == "t"][v$nodehosts[v$nodetypes == "t"]]))
#   newmeanG <- with(phybreak.object,
#                    p$shape.gen / rgamma(1,shape = p$shape.gen * p$obs + prior.sh,
#                                         rate = sumgt + prior.sc/p$shape.gen))
#   return(within(phybreak.object,
#                 p$mean.gen <- newmeanG))
#
# }
#
#
