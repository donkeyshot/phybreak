add_history <- function(d, vars, pars, h, s, build = FALSE){
  if (build){
    v <- vars
    p <- pars
    copy2pbe1_2("d", environment())
    copy2pbe1_2("v", environment())
    copy2pbe1_2("p", environment())
    copy2pbe1_2("h", environment())
  }
  
  introductions <- sum(vars$infectors==0)
  if (introductions == 1 & build){
    buildoneintro <- TRUE
  #   introductions <- 2
  } else {
    buildoneintro <- FALSE
  }

  h$dist <- rbind(0.5, cbind(0.5, h$dist))
  
  nodetimes <- vars$nodetimes
  nodeparents <- as.integer(vars$nodeparents)
  nodehosts <- as.integer(vars$nodehosts) + 1
  nodetypes <- vars$nodetypes
  inftimes <- c(-1e10, vars$inftimes)
  infectors <- c(0, vars$infectors + 1)
  nhosts <- length(inftimes)
  nsamples <- sum(nodetypes %in% c("s", "x"))
  ncoalnodes <- sum(nodetypes == "c")
  ntransnodes <- sum(nodetypes == "t")
  
  names(nodetimes) <- c(rep("s", nsamples), 
                        rep("c", ncoalnodes),
                        rep("t", ntransnodes))
  names(nodeparents) <- c(rep("s", nsamples), 
                        rep("c", ncoalnodes),
                        rep("t", ntransnodes))
  
  # artificial host
  nodehosts <- c(nodehosts[1:nsamples], 
                 rep(1, introductions -1), nodehosts[nodetypes == "c"],
                 0, nodehosts[nodetypes == "t"])
  
  # artificial types
  nodetypes <- c(rep("s", nsamples),
                 rep("c", nsamples - 1),
                 rep("t", ntransnodes + 1))
    
  
  # artificial times
  edgesin <- which(nodehosts == 1 & nodetypes != "c")
  edgeintimes <- inftimes[infectors == 1]
  # if (oneintro) {
  #   edgesin <- rep(edgesin, 2)
  #   edgeintimes <- rep(edgeintimes, 2)
  # }
  newcoaltimes <- sample_coaltimes(edgeintimes, inftimes[1], pars)
  
  nodetimes <- c(nodetimes[1:nsamples], 
                 newcoaltimes, 
                 nodetimes[names(nodetimes) == "c"],
                 inftimes[1],
                 nodetimes[names(nodetimes) == "t"])

  # artificial parents
  coalnodes <- which(nodetypes == "c")[1:length(newcoaltimes)]
  transnode <- which(nodetypes == "t")[1]
  nodeorder <- order(c(newcoaltimes, edgeintimes))
  edgeend <- c(coalnodes, edgesin)[nodeorder]
  edgeendtimes <- c(newcoaltimes, edgeintimes)[nodeorder]
  
  # sample topology of minitree within hostID
  edgestart <- sample_topology(edgeend,
                               edgeendtimes,
                               c(rep("c", length(newcoaltimes)),
                                 rep("x", length(edgeintimes)))[nodeorder],
                               transnode)
  
  
  nodeparents <- sapply(nodeparents, function(xx){
    if (xx == 0) xx <- 1
    else if (xx <= sum(nsamples + ncoalnodes)) xx <- xx + length(newcoaltimes)
    else xx <- xx + length(newcoaltimes) + 1
  })
  nodeparents <- c(nodeparents[names(nodeparents) == "s"],
                        rep(NA, length(newcoaltimes)),
                        nodeparents[names(nodeparents) == "c"],
                        0,
                        nodeparents[names(nodeparents) == "t"]) 
  
  if (buildoneintro) nodeparents[transnode+which(infectors == 1) -1] <- transnode
  else nodeparents[edgeend] <- edgestart
  # make element v
  v <- list()
  v$inftimes <- inftimes
  v$infectors <- infectors
  v$nodetypes <- nodetypes
  v$nodetimes <- nodetimes
  v$nodehosts <- nodehosts
  v$nodeparents <- nodeparents
  
  if (sum(v$infectors==1) == sum(pbe0$hs$infectors==1)) {
    if (all(which(v$infectors==1) == which(pbe0$hs$infectors==1))) {
      if (all(v$inftimes[v$infectors==1] == pbe0$hs$inftimes[v$infectors==1])) {
        v$inftimes[1] <- pbe0$hs$inftimes[1]
        v$nodetimes[pars$obs+1:length(newcoaltimes)] <- pbe0$hs$nodetimes[pars$obs+1:length(newcoaltimes)]
        v$nodetimes[sum(v$nodetypes!="t")+1] <- v$inftimes[1]
        
        for(n in names(pbe1))
          copy2pbe1_2(n, pbe1)
        
        hs <- v
        copy2pbe1("v", environment())
        copy2pbe1("hs", environment())
        return()
      }
    }
  }
  
  return(list("v" = v, "h" = h))
  
  ctree <- list(d = d, v = environment2phybreak(v), p = pars, h = h, s = s, history = TRUE)
  class(ctree) <- c("phybreak", "list")
  
  for(n in names(pbe0))
    copy2pbe0_2(n, pbe0)
  for(n in names(pbe1))
    copy2pbe1_2(n, pbe1)

  rm(list = ls(pbe0), envir = pbe0)
  rm(list = ls(pbe1), envir = pbe1)
  
  # set the number of cycles
  nsample = 500
  
  build_pbe(ctree)
  print("ctree built")
  copy2pbe1("logLiksam", pbe0)
  copy2pbe1("history", pbe0)
  
  ### create room in s to add the new posterior samples
  s.post <- list(inftimes = with(ctree, cbind(s$inftimes, matrix(NA, nrow = 1, ncol = nsample))),
                 nodetimes = with(ctree, cbind(s$nodetimes, matrix(NA, nrow = length(newcoaltimes), ncol = nsample))),
                 nodeparents = with(ctree, cbind(s$nodeparents, matrix(NA, nrow = length(v$nodeparents)+ p$obs + 1, ncol = nsample))),
                 logLik = c(ctree$s$logLik, rep(NA, nsample)))
  
  # s.post <- list(inftimes = with(x, cbind(s$inftimes, matrix(NA, nrow = p$obs+1, ncol = nsample))),
  #                infectors = with(x, cbind(s$infectors, matrix(NA, nrow = p$obs+1, ncol = nsample))),
  #                nodetimes = with(x, cbind(s$nodetimes, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
  #                nodehosts = with(x, cbind(s$nodehosts, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
  #                nodeparents = with(x, cbind(s$nodeparents, matrix(NA, nrow = 2 * d$nsamples - 1, ncol = nsample))),
  #                introductions = c(sum(x$s$infectors==0), rep(NA, nsample)),
  #                mu = c(x$s$mu, rep(NA, nsample)), 
  #                mG = c(x$s$mG, rep(NA, nsample)), 
  #                mS = c(x$s$mS, rep(NA, nsample)), 
  #                wh.s = c(x$s$wh.s, rep(NA, nsample)), 
  #                wh.e = c(x$s$wh.e, rep(NA, nsample)), 
  #                wh.0 = c(x$s$wh.0, rep(NA, nsample)), 
  #                dist.e = c(x$s$dist.e, rep(NA, nsample)), 
  #                dist.s = c(x$s$dist.s, rep(NA, nsample)), 
  #                dist.m = c(x$s$dist.m, rep(NA, nsample)), 
  #                logLik = c(x$s$logLik, rep(NA, nsample)))
  
  protocoldistribution <- c(1/3, 1/3, 1/3)
  protocoldistribution <- sapply(protocoldistribution, round, digits = 5)
  
  for (sa in tail(1:length(s.post$logLik), nsample)){
    # print(sprintf("sa = %s",sa))
    for (i in sample(c(-1,1))){
      print(i)
      if (i == 1){
        which_protocol <- sample(c("classic", "keepphylo", "withinhost"),
                                 1,
                                 prob = protocoldistribution)
        # print(which_protocol)
        update_host(1, "keepphylo")
        # print(pbe0$v)
      }
      
      if (i == -1 && ctree$h$est.mG)  update_mGh()
    }
    # print(pbe0$p)
    
    s.post$inftimes[, sa] <- pbe0$v$inftimes[1]
    s.post$nodetimes[, sa] <- pbe0$v$nodetimes[pbe0$v$nodetypes == "c"][1:length(newcoaltimes)]
    s.post$nodeparents[, sa] <- pbe0$v$nodeparents
    s.post$logLik[sa] <- pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikdist + pbe0$logLikcoal
    
    # vars_to_log <- environment2phybreak(pbe0$v)
    # s.post$inftimes[, sa] <- vars_to_log$inftimes
    # s.post$infectors[, sa] <- vars_to_log$infectors
    # s.post$nodetimes[, sa] <- c(vars_to_log$nodetimes[vars_to_log$nodetypes == "c"],
    #                             rep(NA,x$d$nsamples-1-length(which(vars_to_log$nodetypes == "c"))))
    # s.post$nodehosts[, sa] <- c(vars_to_log$nodehosts[vars_to_log$nodetypes == "c"],
    #                             rep(NA,x$d$nsamples-1-length(which(vars_to_log$nodetypes == "c"))))
    # s.post$nodeparents[, sa] <- c(vars_to_log$nodeparents,
    #                               rep(NA,x$d$nsamples-1-length(which(vars_to_log$nodetypes == "c"))))
    # s.post$introductions[sa] <- sum(vars_to_log$infectors == 0)
    # s.post$mu[sa] <- pbe0$p$mu
    # s.post$mG[sa] <- pbe0$p$gen.mean
    # s.post$mS[sa] <- pbe0$p$sample.mean
    # s.post$wh.s[sa] <- pbe0$p$wh.slope
    # s.post$wh.e[sa] <- pbe0$p$wh.exponent
    # s.post$wh.0[sa] <- pbe0$p$wh.level
    # s.post$dist.e[sa] <- pbe0$p$dist.exponent
    # s.post$dist.s[sa] <- pbe0$p$dist.scale
    # s.post$dist.m[sa] <- pbe0$p$dist.mean
    # s.post$logLik[sa] <- pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikdist + pbe0$logLikcoal
  }
  
  index <- which(pbe0$v$infectors==1)
  index_trans <- which(pbe0$v$nodetypes == "t")[index]
  np <- apply(s.post$nodeparents, 2, function(xx) paste(xx, collapse = ""))
  np <- names(table(np))[which.max(table(np))]
  np_index <- which(apply(s.post$nodeparents, 2, function(xx) paste(xx, collapse = "")) == np)
  tree <- np_index[which.max(s.post$logLik[np_index])]
  
  # med_logLik <- which.min(abs(median(s.post$logLik[(nsample-500):nsample]) - s.post$logLik[(nsample-500):nsample]))[1]
  # med_logLik <- nsample - 500 + med_logLik
  # 
  hs <- pbe0$v
  hs$inftimes[1] <- s.post$inftimes[tree]
  hs$nodetimes[(pbe0$p$obs + 1):(pbe0$p$obs + length(newcoaltimes))] <- s.post$nodetimes[,tree]
  hs$nodetimes[sum(hs$nodetypes!="t")+1] <- hs$inftimes[1]
  hs$nodeparents <- s.post$nodeparents[,tree]
  
  history <- FALSE
  
  if (build){
    copy2pbe0("history", environment())
    copy2pbe0("hs", environment())
    copy2pbe0("v", pbe1_2)
    
    
  } else {
    rm(list = ls(pbe0), envir = pbe0)
    for(n in names(pbe0_2))
      copy2pbe0(n, pbe0_2)
    
    copy2pbe1("hs", environment())
  }
}

remove_history <- function(keepenv = FALSE){
  v <- pbe0$v
  p <- pbe0$p
  d <- pbe0$d
  
  v$inftimes <- tail(v$inftimes, p$obs)
  v$infectors <- tail(v$infectors, p$obs) - 1
  intros <- sum(v$infectors == 0)
  
  keephosts <- setdiff(which(v$nodehosts != 0), 
                       intersect(which(v$nodehosts == 1), which(v$nodetypes == "c")))
  v$nodehosts <- v$nodehosts[keephosts] - 1
  v$nodetypes <- v$nodetypes[keephosts]
  v$nodetimes <- v$nodetimes[keephosts]
  
  remnodes <- setdiff(1:length(v$nodeparents), keephosts)
  transnode <- 2*d$nsamples - length(remnodes) + 1
  nodeparents <- v$nodeparents
  nodeparents[nodeparents %in% remnodes] <- 0
  for (i in seq_along(remnodes)){
    remnode <- remnodes[i] - i + 1
    nodeparents <- sapply(nodeparents, function(xx) {
      if (xx > remnode) xx - 1 
      else xx
    })
    nodeparents <- nodeparents[-remnode]
  }
  v$nodeparents <- nodeparents
  
  if (keepenv)
    copy2pbe0_2("v", environment())
  else {
    pbe0$h$dist <- pbe0$h$dist[,2:ncol(pbe0$h$dist)][2:nrow(pbe0$h$dist),]
    copy2pbe0("v", environment())
  }
}