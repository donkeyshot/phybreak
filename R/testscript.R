library(phybreak)
simul <- sim_phybreak(10, samplesperhost = 2)
res <- .obsarray(matrix(unlist(simul$sequences), ncol = 20), 20)
dim(res) <- c(15,16,20)
res[,,10]
sum(res)
reshost <- obshostarrayC(res, rep(1:10,2), 10)
dim(reshost) <- c(15,16,10)

tramatR <- make_transitionmatrix(
  attr(simul$sequences,"contrast")[c(1,2,6,3,7,9,12,4,8,10,13,11,14,15,16),],.001,.1)
tramatC <- transitmatrixC(.001,.1)
dim(tramatC) <- c(15,15)

simul <- sim_phybreak(10, samplesperhost = 2)
pbo <- phybreak(simul, use.tree = T, wh.model = "no_coalescent", mu = .001)
likgenetic(matrix(unlist(simul$sequences), ncol = 20), 
           attr(simul$sequences, "weight"), 
           match(simul$sim.infectors, c("index", names(simul$sim.infectors)))-1, 
           match(simul$sample.hosts, c(names(simul$sim.infectors))),
           .001, .1)
