# The environments .phybreakenv and .phybreakenv.prop will contain phybreak objects (excluding mcmc samples) + likelihood
# calculations.  .phybreakenv is consistent with the current state of updating phybreak.object; .phybreakenv.prop is used to
# propose states and calculate likelihoods for these proposals.

# likarray is an array with likelihoods per nucleotide per SNP per node (sampling and coalescence). They are stored as
# vectors of length Nnodes*nSNPs*4, with implicit dim = c(4,nSNPs,Nnodes).

# logLikgen and logLiksam are the log-likelihood values of generation and sampling times with means mG and mS. logLikcoal is
# the log-likelihood value of the coalescent model.

# The function '.build.phybreakenv' is used to initialise The function '.destroy.phybreakenv' is used at the end The function
# '.prepare.phybreakenv' is used to prepare .phybreakenv.prop for another proposal The function '.propose.phybreakenv' is
# used to calculate the likelihoods in .phybreakenv.prop after changing in proposal functions.  The function
# '.accept.phybreakenv' is to change .phybreakenv by copying .phybreakenv.prop.


# The environments are only used during MCMC-updating.
.phybreakenv <- new.env()
.phybreakenv.prop <- new.env()


### build the environments at the start of an mcmc chain by copying the fixed parameters and phybreak object, and by
### calculating likarray and the log-likelihoods called from: burnin.phybreak sample.phybreak calls: .likseqenv ## C++ function
### .liksampletimes .likgentimes .likcoaltimes .accept.phybreakenv
.build.phybreakenv <- function(phybreak.obj) {
    .phybreakenv$tinf.prop.shape.mult <- tinf.prop.shape.mult
    .phybreakenv.prop$tinf.prop.shape.mult <- tinf.prop.shape.mult
    
    .phybreakenv.prop$d <- phybreak.obj$d
    .phybreakenv.prop$v <- phybreak.obj$v
    .phybreakenv.prop$p <- phybreak.obj$p
    .phybreakenv.prop$h <- phybreak.obj$h
    
    ### Jukes-Cantor: reduce diversity by naming each nucleotide by its frequency order, and grouping SNPs by same pattern across
    ### hosts
    fn <- function(snpvector) {
        match(snpvector, names(sort(table(snpvector), decreasing = TRUE)))
    }
    snpreduced <- apply(phybreak.obj$d$SNP, 2, fn)
    snpfrreduced <- phybreak.obj$d$SNPfr
    snpreduced[phybreak.obj$d$SNP == "n"] <- 5
    if (ncol(phybreak.obj$d$SNP) > 1) {
        for (i in (ncol(phybreak.obj$d$SNP) - 1):1) {
            for (j in length(snpfrreduced):(i + 1)) {
                if (all(snpreduced[, i] == snpreduced[, j])) {
                  snpfrreduced[i] <- snpfrreduced[i] + snpfrreduced[j]
                  snpfrreduced <- snpfrreduced[-j]
                  snpreduced <- snpreduced[, -j, drop = FALSE]
                }
            }
        }
    }
    .phybreakenv.prop$likarrayfreq <- snpfrreduced
    
    
    ### initialize all dimensions of likarray
    .phybreakenv.prop$likarray <- array(1, dim = c(4, length(snpfrreduced), 2 * .phybreakenv.prop$p$obs - 1))
    ### initialize likarray with observations on sampling nodes: 0 or 1
    .phybreakenv.prop$likarray[cbind(1, rep(1:length(snpfrreduced), .phybreakenv.prop$p$obs), rep(1:.phybreakenv.prop$p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 1 | snpreduced == 5)
    .phybreakenv.prop$likarray[cbind(2, rep(1:length(snpfrreduced), .phybreakenv.prop$p$obs), rep(1:.phybreakenv.prop$p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 2 | snpreduced == 5)
    .phybreakenv.prop$likarray[cbind(3, rep(1:length(snpfrreduced), .phybreakenv.prop$p$obs), rep(1:.phybreakenv.prop$p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 3 | snpreduced == 5)
    .phybreakenv.prop$likarray[cbind(4, rep(1:length(snpfrreduced), .phybreakenv.prop$p$obs), rep(1:.phybreakenv.prop$p$obs, 
        each = length(snpfrreduced)))] <- 1 * t(snpreduced == 4 | snpreduced == 5)
    
    
    ### complete likarray and calculate log-likelihood of sequences
    .likseqenv(.phybreakenv.prop, (.phybreakenv.prop$p$obs + 1):(2 * .phybreakenv.prop$p$obs - 1), 1:.phybreakenv.prop$p$obs)
    
    ### calculate the other log-likelihoods
    .phybreakenv.prop$logLiksam <- with(phybreak.obj, .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes, v$nodetypes))
    
    .phybreakenv.prop$logLikgen <- with(phybreak.obj, .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes, v$nodehosts, 
        v$nodetypes))
    
    .phybreakenv.prop$logLikcoal <- with(phybreak.obj, .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, 
        v$nodetypes))
    
    ### copy everything into .phybreakenv
    .accept.phybreakenv()
}


### take the elements d, v, p, and h from .phybreakenv, and s from the function arguments, and make a new phybreak-object. Then
### empty the environments and return the new object.  called from: burnin.phybreak sample.phybreak
.destroy.phybreakenv <- function(phybreak.obj.samples) {
    res <- list(d = .phybreakenv$d, v = .phybreakenv$v, p = .phybreakenv$p, h = .phybreakenv$h, s = phybreak.obj.samples)
    class(res) <- c("phybreak", "list")
    rm(list = ls(.phybreakenv), envir = .phybreakenv)
    rm(list = ls(.phybreakenv.prop), envir = .phybreakenv.prop)
    return(res)
}


### copy the elements from .phybreakenv into .phybreakenv.prop to prepare for a proposal called from: .updatehost .update.mu,
### .update.mS, .update.mG, .update.wh
.prepare.phybreakenv <- function() {
    .phybreakenv.prop$d <- .phybreakenv$d
    .phybreakenv.prop$v <- .phybreakenv$v
    .phybreakenv.prop$p <- .phybreakenv$p
    .phybreakenv.prop$h <- .phybreakenv$h
    .phybreakenv.prop$likarray <- .phybreakenv$likarray + 0  #make a true copy, not a pointer
    .phybreakenv.prop$likarrayfreq <- .phybreakenv$likarrayfreq
    .phybreakenv.prop$logLik <- .phybreakenv$logLik
    .phybreakenv.prop$logLiksam <- .phybreakenv$logLiksam
    .phybreakenv.prop$logLikgen <- .phybreakenv$logLikgen
    .phybreakenv.prop$logLikcoal <- .phybreakenv$logLikcoal
}


### calculate the new log-likelihoods where necessary and adjust likarray. Argument f indicates which type of function it is
### called from called from: .updatepathA - .updatepathJ .update.mu, .update.mS, .update.mG, .update.wh calls: .likseqenv ##
### C++ function .liksampletimes .likgentimes .likcoaltimes
.propose.phybreakenv <- function(f) {
    if (f == "phylotrans") {
        # identify changed nodes
        chnodes <- which((.phybreakenv.prop$v$nodeparents != .phybreakenv$v$nodeparents) | (.phybreakenv.prop$v$nodetimes != 
            .phybreakenv$v$nodetimes))
        chnodes <- unique(unlist(sapply(chnodes, .ptr, pars = .phybreakenv.prop$v$nodeparents)))
        chnodes <- chnodes[chnodes > .phybreakenv.prop$p$obs & chnodes < 2 * .phybreakenv.prop$p$obs]
        # identify nodetips
        nodetips <- c(match(chnodes, .phybreakenv.prop$v$nodeparents), 3 * .phybreakenv.prop$p$obs - match(chnodes, rev(.phybreakenv.prop$v$nodeparents)))
        nodetips[nodetips >= 2 * .phybreakenv.prop$p$obs] <- match(nodetips[nodetips >= 2 * .phybreakenv.prop$p$obs], .phybreakenv.prop$v$nodeparents)
        nodetips <- nodetips[is.na(match(nodetips, chnodes))]
    } else if (f == "mu") {
        chnodes <- (.phybreakenv.prop$p$obs + 1):(2 * .phybreakenv.prop$p$obs - 1)
        nodetips <- 1:.phybreakenv.prop$p$obs
    } else {
        chnodes <- NULL
    }
    
    
    if (!is.null(chnodes)) {
        .likseqenv(.phybreakenv.prop, chnodes, nodetips)
    }
    
    
    if (f == "phylotrans" || f == "trans" || f == "mG") {
        evalq({
            logLikgen <- .lik.gentimes(p$obs, p$shape.gen, p$mean.gen, v$nodetimes, v$nodehosts, v$nodetypes)
        }, envir = .phybreakenv.prop)
    }
    
    if (f == "phylotrans" || f == "trans" || f == "mS") {
        evalq({
            logLiksam <- .lik.sampletimes(p$shape.sample, p$mean.sample, v$nodetimes, v$nodetypes)
        }, envir = .phybreakenv.prop)
    }
    
    if (f == "phylotrans" || f == "trans" || f == "slope") {
        evalq({
            logLikcoal <- .lik.coaltimes(p$obs, p$wh.model, p$wh.slope, v$nodetimes, v$nodehosts, v$nodetypes)
        }, envir = .phybreakenv.prop)
    }
    
}


### copy the elements from .phybreakenv.prop into .phybreakenv upon acceptance of a proposal called from: .updatepathA -
### .updatepathJ .update.mu, .update.mS, .update.mG, .update.wh
.accept.phybreakenv <- function() {
    .phybreakenv$d <- .phybreakenv.prop$d
    .phybreakenv$v <- .phybreakenv.prop$v
    .phybreakenv$p <- .phybreakenv.prop$p
    .phybreakenv$h <- .phybreakenv.prop$h
    .phybreakenv$likarray <- .phybreakenv.prop$likarray
    .phybreakenv$likarrayfreq <- .phybreakenv.prop$likarrayfreq
    .phybreakenv$logLik <- .phybreakenv.prop$logLik
    .phybreakenv$logLiksam <- .phybreakenv.prop$logLiksam
    .phybreakenv$logLikgen <- .phybreakenv.prop$logLikgen
    .phybreakenv$logLikcoal <- .phybreakenv.prop$logLikcoal
}

