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
    v <- phybreak.obj$v
    p <- phybreak.obj$p
    SNP <- t(matrix(unlist(d$sequences), ncol = p$obs))
    SNPfr <- attr(d$sequences, "weight")
  
    
    ### Jukes-Cantor: reduce diversity by naming each nucleotide by its frequency order, and grouping SNPs by same pattern across
    ### hosts
    fn <- function(snpvector) {
        match(snpvector, names(sort(table(snpvector), decreasing = TRUE)))
    }
    snpreduced <- apply(SNP, 2, fn)
    snpfrreduced <- SNPfr
    snpreduced[SNP == 16] <- 5
    if (ncol(SNP) > 1) {
        for (i in (ncol(SNP) - 1):1) {
            for (j in length(snpfrreduced):(i + 1)) {
                if (all(snpreduced[, i] == snpreduced[, j])) {
                  snpfrreduced[i] <- snpfrreduced[i] + snpfrreduced[j]
                  snpfrreduced <- snpfrreduced[-j]
                  snpreduced <- snpreduced[, -j, drop = FALSE]
                }
            }
        }
    }
    likarrayfreq <- snpfrreduced
    
    
    ### initialize all dimensions of likarray
    likarray <- array(1, dim = c(4, length(snpfrreduced), 2 * p$obs - 1))
    ### initialize likarray with observations on sampling nodes: 0 or 1
    likarray[cbind(1, rep(1:length(snpfrreduced), p$obs), rep(1:p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 1 | snpreduced == 5)
    likarray[cbind(2, rep(1:length(snpfrreduced), p$obs), rep(1:p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 2 | snpreduced == 5)
    likarray[cbind(3, rep(1:length(snpfrreduced), p$obs), rep(1:p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 3 | snpreduced == 5)
    likarray[cbind(4, rep(1:length(snpfrreduced), p$obs), rep(1:p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 4 | snpreduced == 5)
    
    
    ### complete likarray and calculate log-likelihood of sequences
    .likseqenv(le, (p$obs + 1):(2 * p$obs - 1), 1:p$obs)

    
        
    ### calculate the other log-likelihoods
    logLiksam <- .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes, v$nodetypes)
    logLikgen <- .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes, v$nodehosts, 
        v$nodetypes)
    logLikcoal <- .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, 
        v$nodetypes)
    
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
    .pbe1$likarray <- .pbe0$likarray + 0  #make a true copy, not a pointer
    .copy2pbe1("likarrayfreq", .pbe0)
}


### calculate the new log-likelihoods where necessary and adjust likarray. Argument f indicates which type of function it is
### called from called from: .updatepathA - .updatepathJ .update.mu, .update.mS, .update.mG, .update.wh calls: .likseqenv ##
### C++ function .liksampletimes .likgentimes .likcoaltimes
.propose.pbe <- function(f) {
    ### Making variables and parameters available within the function
    le <- environment()
    v <- .pbe1$v
    p <- .pbe1$p
  
    if (f == "phylotrans") {
        # identify changed nodes
        chnodes <- which((v$nodeparents != .pbe0$v$nodeparents) | (v$nodetimes != 
            .pbe0$v$nodetimes))
        chnodes <- unique(unlist(sapply(chnodes, .ptr, pars = v$nodeparents)))
        chnodes <- chnodes[chnodes > p$obs & chnodes < 2 * p$obs]
        # identify nodetips
        nodetips <- c(match(chnodes, v$nodeparents), 3 * p$obs - match(chnodes, rev(v$nodeparents)))
        nodetips[nodetips >= 2 * p$obs] <- match(nodetips[nodetips >= 2 * p$obs], v$nodeparents)
        nodetips <- nodetips[is.na(match(nodetips, chnodes))]
    } else if (f == "mu") {
        chnodes <- (p$obs + 1):(2 * p$obs - 1)
        nodetips <- 1:p$obs
    } else {
        chnodes <- NULL
    }
    
    
    if (!is.null(chnodes)) {
        .likseqenv(.pbe1, chnodes, nodetips)
    }
    
    
    if (f == "phylotrans" || f == "trans" || f == "mG") {
      logLikgen <- .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes, v$nodehosts, v$nodetypes)
      .copy2pbe1("logLikgen", le)
    }
    
    if (f == "phylotrans" || f == "trans" || f == "mS") {
      logLiksam <- .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes, v$nodetypes)
      .copy2pbe1("logLiksam", le)
    }
    
    if (f == "trans" || f == "slope") {
      logLikcoal <- .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes)
      .copy2pbe1("logLikcoal", le)
    }
    
}


### copy the elements from .pbe1 into .pbe0 upon acceptance of a proposal called from: .updatepathA -
### .updatepathJ .update.mu, .update.mS, .update.mG, .update.wh
.accept.pbe <- function(f) {
    .copy2pbe0("v", .pbe1)
    .copy2pbe0("p", .pbe1)
    
    if(f == "phylotrans" || f == "mu") {
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
      logLikcoal <- with(.pbe1, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes))
      .copy2pbe0("logLikcoal", environment())
    }
    
    rm(list = ls(.pbe1), envir = .pbe1)
}

