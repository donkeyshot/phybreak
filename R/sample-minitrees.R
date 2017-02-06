### functions to simulate mini-trees ###



### sample coalescent times in a host, given the tip times since infection,
### ...the wh-model and the slope (used if WHmodel = 3)
.samplecoaltimes <- function(tleaves, WHmodel = 3, slope = 1, reassortment = 0, Ngenes = 1) {
  ### tests
  if(min(tleaves) < 0) stop(".samplecoaltimes with negative tip times")
  if(!any(WHmodel == 1:3)) stop(paste0(".samplecoaltimes called with WHmodel = ",WHmodel))
  if(WHmodel == 3 && slope < 0) stop(".samplecoaltimes called with negative slope")

  ### function body
  if(length(tleaves) < 2) return(c()) # May give an error
  
  switch(WHmodel,
         # coalescence at transmission
         return(matrix(rep(head(sort(tleaves),-1),Ngenes), nrow = Ngenes), byrow = TRUE),
         # coalescence at infection
         return(matrix(rep(0*tleaves[-1],Ngenes), nrow = Ngenes), byrow = TRUE),
         { # transform times so that fixed rate 1 can be used
           ttrans <- sort(log(tleaves)/(slope), decreasing = TRUE)
           if (reassortment != 0){
             tnodetrans <- matrix(sapply(1:Ngenes, function(gene) .sctwh3(ttrans)), nrow = Ngenes, byrow = TRUE)
           } else {
             tnodetrans <- matrix(rep(.sctwh3(ttrans), Ngenes), nrow = Ngenes, byrow = TRUE)
           }
           res <- matrix(sapply(1:Ngenes, function(gene) sort(exp(slope*tnodetrans[gene, ]))), nrow = Ngenes, byrow = TRUE)
           
           # make sure that all branches will have positive length
           res <- matrix(
             sapply(1:Ngenes, function(gene) apply(cbind(res[gene,], min(10^-5,tleaves/length(tleaves))
                                                         *(1:length(res[gene,]))),1,max)), nrow = Ngenes, byrow = TRUE)
           return(res)
         }
  )
}


### sample tree topology in a host, given the node IDs, ntimes, and types,
### ...the root node and the WHmodel
.sampletopology <- function(nIDs, ntimes, ntypes, rootnode, WHmodel = 3, reassortment, Ngenes) {
  ### tests
  if(!any(WHmodel == 1:3)) stop(paste0(".sampletopology called with WHmodel = ",WHmodel))
  ### function body
  if(length(nIDs) == 1) return(rootnode)
  
  ntimes <- matrix(ntimes, nrow = Ngenes)  # Required for gene == 1 case 

  switch(
    WHmodel,
    #coalescence at transmission
    {
      cnodes <- nIDs[ntypes=="c"]
      cnodeparents <- c(rootnode,head(cnodes,-1))
      leafparents <- c(cnodes,tail(cnodes,1))
      leafparents <- matrix(sapply(1:Ngenes, function(gene) leafparents[rank(ntimes[gene, ntypes != "c"],ties.method="first")]),
                            nrow = Ngenes, byrow = TRUE)
      res<- matrix( sapply(1:Ngenes, function(gene) 
        c(head(leafparents[Ngenes, ],sum(ntypes=="s")),cnodeparents, tail(leafparents[Ngenes, ],sum(ntypes=="t")))), 
        nrow = Ngenes, byrow = TRUE)
      return(res)
    },
    #coalescence at infection
    {
      cnodes <- sort(nIDs[ntypes=="c"],decreasing = TRUE)
      res <- c(rep(NA,sum(ntypes=="s")),
               rootnode,tail(-nIDs[ntypes=="c"],-1),
               rep(NA,sum(ntypes=="t")))
      for(i in cnodes) {
        res[sample(which(is.na(res)),2)] <- i
        res[res == -i] <- NA
      }
      return(matrix(rep(res, Ngenes), nrow = dim(ntimes)[1], byrow = TRUE))
    }
  )
  
  res2 <- c()
  for(ii in 1:(reassortment*(Ngenes-1)+1)){
    IDs <- nIDs[order(ntimes[ii, ],ntypes)]
    tys <- ntypes[order(ntimes[ii, ],ntypes)]
    if(tys[1] != "c") {
      print(c(nIDs, ntimes[ii,], ntypes, rootnode))
      stop("host topology does not start with coalescence node")
    }
    res <- rep(rootnode, length(nIDs))
    tochoose <- rep(IDs[1], 2)
    for(i in 2:length(nIDs)) {
      res[i] <- tochoose[1]
      if(tys[i] == "c") {
        tochoose <- sample(c(tochoose[-1], IDs[i], IDs[i]))
      } else {
        tochoose <- tochoose[-1]
      }
    }
    res2 <- rbind(res2, res[order(IDs)])
  }
  if(reassortment == FALSE){res2 <- matrix(rep(res2,Ngenes),nrow = Ngenes, byrow = TRUE)}
  
  return(res2)
}