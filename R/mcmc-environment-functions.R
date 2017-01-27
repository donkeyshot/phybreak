
# The environments .pbe0 and .pbe1 will contain phybreak objects (excluding mcmc samples) + likelihood
# calculations.  .pbe0 is consistent with the current state of updating phybreak.object; .pbe1 is used to
# propose states and calculate likelihoods for these proposals.

# likarray is an array with likelihoods per nucleotide per SNP per node (sampling and coalescence). They are stored as
# vectors of length Nnodes*nSNPs*4, with implicit dim = c(4,nSNPs,Nnodes).

# logLikgen and logLiksam are the log-likelihood values of generation and sampling times with means mG and mS. logLikcoal is
# the log-likelihood value of the coalescent model.

# The function '.build.pbe' is used to initialise The function '.destroy.pbe' is used at the end The function
# '.prepare.pbe' is used to prepare .pbe1 for another proposal The function '.propose.pbe' is
# used to calculate the likelihoods in .pbe1 after changing in proposal functions.  The function
# '.accept.pbe' is to change .pbe0 by copying .pbe1.


# The environments are only used during MCMC-updating.
.pbe0 <- new.env()
.pbe1 <- new.env()

# Copy functions to phybreak environments
.copy2pbe0 <- function(var, env) {
  assign(var, get(var, env), .pbe0)
}
.copy2pbe1 <- function(var, env) {
  assign(var, get(var, env), .pbe1)
}


### build the .pbe0 at the start of an mcmc chain by copying the fixed parameters and phybreak object, and by
### calculating likarray and the log-likelihoods called from: burnin.phybreak sample.phybreak calls: .likseqenv ## C++ function
### .liksampletimes .likgentimes .likcoaltimes .accept.pbe
.build.pbe <- function(phybreak.obj) {
  ### Making everything available within the function
  
  le <- environment()
  d <- phybreak.obj$d
  h <- phybreak.obj$h
  v <- phybreak.obj$v  # Later defined again in the if branch, fix this.
  p <- phybreak.obj$p
  Ngenes = length(d$nSNPs)
  SNP <- lapply(1:Ngenes, function(x) t(matrix(unlist(d$sequences[[x]]), ncol =p$obs)))
  SNPfr <- lapply(1:Ngenes,function(x) attr(d$sequences[[x]], "weight"))              
  GeneNames <- paste("gene",1:Ngenes,sep="");    names(SNPfr)<- GeneNames; names(SNP) <- GeneNames # Doe ik hier nog wat mee, zo niet WEG!
  ### Jukes-Cantor: reduce diversity by naming each nucleotide by its frequency order, and grouping SNPs by same pattern across
  ### hosts
  
  fn <- function(snpvector) {match(snpvector, names(sort(table(snpvector), decreasing = TRUE)))}
  snpreduced <- lapply(1:Ngenes, function(x) {apply(SNP[[x]], 2, fn)}); names(snpreduced) = GeneNames
  snpfrreduced <- SNPfr
  for(i in 1:Ngenes){ snpreduced[[i]][  SNP[[i]] == 16  ] <- 5 }  
  
  for(gene in 1:Ngenes){
    if (ncol(SNP[[gene]]) > 1) {
      for (i in (ncol(SNP[[gene]]) - 1):1) {
        for (j in length(snpfrreduced[[gene]]):(i + 1)) {
          if (all(snpreduced[[gene]][, i] == snpreduced[[gene]][, j])) {
            snpfrreduced[[gene]][i] <- snpfrreduced[[gene]][i] + snpfrreduced[[gene]][j]
            snpfrreduced[[gene]] <- snpfrreduced[[gene]][-j]
            snpreduced[[gene]] <- snpreduced[[gene]][, -j, drop = FALSE]
          }
        }
      }
    }
  }
  likarrayfreq    <- snpfrreduced  
  Templikarrayfreq <- likarrayfreq
  
  # Initialize collection list and vectors
  LikCollect <- list() 
  LikCollect <- sapply(paste("gene", 1:Ngenes, sep = ""), function(x) NULL)
  logLikcoal <- rep(0,Ngenes)
  logLikseqCollect <- rep(0,Ngenes)
  
  for(gene in 1:Ngenes){
    ### initialize all dimensions of likarray
    likarray <- array(1, dim = c(4, length(snpfrreduced[[gene]]), 2 * p$obs - 1)) # 3D array
    likarrayfreq <- Templikarrayfreq[[gene]]
    
    # Note that nodehost will probably become a matrix.
    v <- list(nodetimes = phybreak.obj$v$nodetimes[gene, ], nodehosts = phybreak.obj$v$nodehosts,
              nodeparents = phybreak.obj$v$nodeparents[gene, ], nodetypes = phybreak.obj$v$nodetypes)
    
    
    ### initialize likarray with observations on sampling nodes: 0 or 1
    likarray[cbind(1, rep(1:length(snpfrreduced[[gene]]), p$obs), rep(1:p$obs, 
                                                                      each = length(snpfrreduced[[gene]])))] <- 1 * t(snpreduced[[gene]] == 1 | snpreduced[[gene]] == 5)
    likarray[cbind(2, rep(1:length(snpfrreduced[[gene]]), p$obs), rep(1:p$obs, 
                                                                      each = length(snpfrreduced[[gene]])))] <- 1 * t(snpreduced[[gene]] == 2 | snpreduced[[gene]] == 5)
    likarray[cbind(3, rep(1:length(snpfrreduced[[gene]]), p$obs), rep(1:p$obs, 
                                                                      each = length(snpfrreduced[[gene]])))] <- 1 * t(snpreduced[[gene]] == 3 | snpreduced[[gene]] == 5)
    likarray[cbind(4, rep(1:length(snpfrreduced[[gene]]), p$obs), rep(1:p$obs, 
                                                                      each = length(snpfrreduced[[gene]])))] <- 1 * t(snpreduced[[gene]] == 4 | snpreduced[[gene]] == 5)
    
    ### complete likarray and calculate log-likelihood of sequences
    .likseqenv(le, (p$obs + 1):(2 * p$obs - 1), 1:p$obs)
    
    
    
    LikCollect[[gene]] <- likarray 
    logLikseqCollect[gene] <- logLikseq
    Templikarrayfreq[[gene]] <- likarrayfreq
  }
  v <- phybreak.obj$v
  
  ### calculate the other log-likelihoods
  logLiksam <- .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes, v$nodetypes) 
  logLikgen <- .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes, v$nodehosts, v$nodetypes) 
  logLikcoal<- .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes, v$reassortment) 
  
  logLikseq <- logLikseqCollect
  likarray <- LikCollect
  likarrayfreq <- Templikarrayfreq
  
  ### copy everything into .pbe0
  .copy2pbe0("d", le)
  .copy2pbe0("h", le)
  .copy2pbe0("v", le)
  .copy2pbe0("p", le)
  .copy2pbe0("likarrayfreq", le)
  .copy2pbe0("likarray", le)
  .copy2pbe0("logLikseq", le)
  .copy2pbe0("logLiksam", le)
  .copy2pbe0("logLikgen", le)
  .copy2pbe0("logLikcoal", le)
}

### take the elements d, v, p, and h from .pbe0, and s from the function arguments, and make a new phybreak-object. Then
### empty the environments and return the new object.  called from: burnin.phybreak sample.phybreak
.destroy.pbe <- function(phybreak.obj.samples) {
  res <- list(d = .pbe0$d, v = .pbe0$v, p = .pbe0$p, h = .pbe0$h, s = phybreak.obj.samples)
  class(res) <- c("phybreak", "list")
  rm(list = ls(.pbe0), envir = .pbe0)
  rm(list = ls(.pbe1), envir = .pbe1)
  return(res)
}

### copy the elements from .pbe0 into .pbe1 to prepare for a proposal called from: .updatehost .update.mu,
### .update.mS, .update.mG, .update.wh
.prepare.pbe <- function() {
  .copy2pbe1("d", .pbe0)
  .copy2pbe1("v", .pbe0)
  .copy2pbe1("p", .pbe0)
  
  .pbe1$likarray <- lapply(1:length(.pbe0$d$nSNPs), function(gene) .pbe0$likarray[[gene]] +0) # make a true copy, not a pointer
  .copy2pbe1("likarrayfreq", .pbe0)
  .pbe1$logLikseq <- lapply(1:length(.pbe0$d$nSNPs), function(gene).pbe0$logLikseq[[gene]] + 0) 
}

## calculate the new log-likelihoods where necessary and adjust likarray. Argument f indicates which type of function it is
### called from called from: .updatepathA - .updatepathJ .update.mu, .update.mS, .update.mG, .update.wh calls: .likseqenv ##
### C++ function .liksampletimes .likgentimes .likcoaltimes
.propose.pbe <- function(f) {
  le <- environment()
  v <- .pbe1$v
  p <- .pbe1$p
  hostID <- .pbe0$hostID
  Ngenes <- dim(v$nodetimes)[1]
  chnodes <- vector(mode = "list", length = Ngenes)
  nodetips <- chnodes
  
  # If reassortment == 1, different number of chnodes may be chosen per gene and
  # nodetips order may vary between genes and sometimes different nodes are chosen!
  if (f == "phylotrans" || f == "topology") {
    for (gene in 1:Ngenes) {
      # identify changed nodes              #
      chnodes[[gene]] <- which((v$nodeparents[gene, ] != .pbe0$v$nodeparents[gene, ]) |
                                 (v$nodetimes[gene, ] != .pbe0$v$nodetimes[gene, ]))
      chnodes[[gene]] <- unique(unlist(sapply(chnodes[[gene]], .ptr, pars = v$nodeparents[gene, ])))
      chnodes[[gene]] <- chnodes[[gene]][chnodes[[gene]] > p$obs & chnodes[[gene]] < 2 * p$obs]
      
      # identify nodetips
      nodetips[[gene]] <- c(match(chnodes[[gene]], v$nodeparents[gene, ]), 
                            3 * p$obs - match(chnodes[[gene]], rev(v$nodeparents[gene, ])))
      nodetips[[gene]][nodetips[[gene]] >= 2 * p$obs] <- match(nodetips[[gene]][nodetips[[gene]] >= 2 * p$obs],
                                                               v$nodeparents[gene, ])
      nodetips[[gene]] <- nodetips[[gene]][is.na(match(nodetips[[gene]], chnodes[[gene]]))]
    }
   # if( (sum(unlist(lapply(nodetips,sum))) == sum(nodetips[[1]])*3) == FALSE){
  #    print("nodetips: FALSE")
  #    print(sum(unlist(lapply(nodetips,sum))) == sum(nodetips[[1]])*3 )
  #    print(nodetips)
    #}
    #if( (sum(unlist(lapply(chnodes,sum))) == sum(chnodes[[1]])*3) == FALSE){
    #  print("chnodes: FALSE")
    #  print(chnodes)
    #}

  } else if (f == "mu") {
    for(gene in 1:Ngenes) {chnodes[[gene]] <-(p$obs + 1):(2 * p$obs - 1)}
    for(gene in 1:Ngenes) {nodetips[[gene]] <- 1:p$obs}
  } else {
    chnodes <- vector(mode = "list", Ngenes)
  }
  
  # Initialize collection list and vectors
  LikCollect <- list()
  logLikseqCollect <- rep(0,Ngenes)
  for (gene in 1:Ngenes) {
    if (!is.null(chnodes[[gene]])) {
      v <- list(nodetimes = .pbe1$v$nodetimes[gene, ], nodehosts = .pbe1$v$nodehosts,
                nodeparents = .pbe1$v$nodeparents[gene, ] , nodetypes = .pbe1$v$nodetypes)
      d <- list(names = .pbe1$d$names, sequences = .pbe1$d$sequences[[gene]],
                sample.times = .pbe1$d$sample.times, nSNPs =.pbe1$d$nSNPS[gene] )
      
      likarray <- .pbe1$likarray[[gene]]
      likarrayfreq <- .pbe1$likarrayfreq[[gene]]
      
      .likseqenv(le, chnodes[[gene]], nodetips[[gene]])
      
      LikCollect[[gene]] <- likarray 
      logLikseqCollect[gene] <- logLikseq
    }
  }
  #print(logLikseqCollect)
  likarray <- LikCollect
  logLikseq <- logLikseqCollect
  
  #rm(LikCollect, logLikseqCollect, envir = le) #, nodeparents, nodetimes, envir = le) 
  .copy2pbe1("logLikseq", le)
  .copy2pbe1("likarray", le)  
  v <- .pbe1$v
  d <- .pbe1$d
  
  if (f == "phylotrans" || f == "trans" || f == "mG") {
    logLikgen <- .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes[1, ], v$nodehosts, v$nodetypes)
    .copy2pbe1("logLikgen", le)
  }
  
  if (f == "phylotrans" || f == "trans" || f == "mS") {
    logLiksam <- .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes[1, ], v$nodetypes)
    .copy2pbe1("logLiksam", le)
  }
  
  if (f == "trans" || f == "slope") {
    logLikcoal <- with(.pbe1, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes,v$reassortment))
    .copy2pbe1("logLikcoal", le)
  }
}


### copy the elements from .pbe1 into .pbe0 upon acceptance of a proposal called from: .updatepathA -
### .updatepathJ .update.mu, .update.mS, .update.mG, .update.wh
.accept.pbe <- function(f) {
  .copy2pbe0("v", .pbe1)
  .copy2pbe0("p", .pbe1)
  
  if(f == "phylotrans" || f == "topology" || f == "mu") {
    .copy2pbe0("likarray", .pbe1)
    .copy2pbe0("logLikseq", .pbe1)
  }
  
  if(f == "phylotrans" || f == "trans" || f == "mS") {
    .copy2pbe0("logLiksam", .pbe1)
  }
  
  if(f == "phylotrans" || f == "trans" || f == "mG") {
    .copy2pbe0("logLikgen", .pbe1)
  }
  
  if(f == "trans" || f == "slope") {
    .copy2pbe0("logLikcoal", .pbe1)
  }
  
  if(f == "phylotrans") {
    logLikcoal <- with(.pbe1, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes,v$reassortment))
    .copy2pbe0("logLikcoal", environment())
  }
  rm(list = ls(.pbe1), envir = .pbe1)
}

