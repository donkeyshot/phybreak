# library(phybreak)
# simul <- sim_phybreak(20, samplesperhost = 2)
# res <- obsarrayC(matrix(unlist(simul$sequences), ncol = 10), 10)
# dim(res) <- c(15,4,10)
# res[,,10]
# sum(res)
# reshost <- obshostarrayC(res, rep(1:5,2), 5)
# dim(reshost) <- c(15,4,5)
# res[,,5]
# 
# tramatR <- make_transitionmatrix(
#   attr(simul$sequences,"contrast")[c(1,2,6,3,7,9,12,4,8,10,13,11,14,15,16),],.001,.1)
# tramatC <- transitmatrixC(.001,.1)
# dim(tramatC) <- c(15,15)
# 
# simul <- sim_phybreak(10, samplesperhost = 2)
# pbo <- phybreak(simul, use.tree = T, wh.model = "no_coalescent", mu = .001)
# system.time({for(i in 1:100) logLik(pbo,T,F,F,F)})
# system.time({for(i in 1:100) likgenetic(matrix(unlist(simul$sequences), ncol = 40),
#            attr(simul$sequences, "weight"),
#            match(simul$sim.infectors, c("index", names(simul$sim.infectors)))-1,
#            match(simul$sample.hosts, c(names(simul$sim.infectors))),
#            1-exp(-.001),1-exp(-.1))})
# 
# 
# 
# simul <- sim_phybreak(20, samplesperhost = 5, mu = 5e-4)
# pbo_linbnw <- phybreak(simul, wh.bottleneck = "wide")
# pbo_linbnw_bi <- burnin_phybreak(pbo_linbnw, 1000)
# pbo_linbnw_sa <- sample_phybreak(pbo_linbnw_bi, 2000)
# ESS(pbo_linbnw_sa)
# transtree(pbo_linbnw_sa, "edmonds")$infector
# 
# pbo_newbnw <- phybreak(simul, wh.model = "no_coalescent", wh.bottleneck = "wide")
# pbo_newbnw_bi <- burnin_phybreak(pbo_newbnw, 1000)
# pbo_newbnw_sa <- sample_phybreak(pbo_newbnw_bi, 2000)
# ESS(pbo_newbnw_sa)
# transtree(pbo_newbnw_sa, "edmonds")$infector
