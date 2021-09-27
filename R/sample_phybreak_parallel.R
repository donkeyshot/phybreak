#' Sampling from a phybreak MCMC-chain.
#' 
#' Function to take (additional) samples from the posterior distribution of a phylogenetic and transmission tree 
#'   (plus associated parameters), within a \code{phybreak} object (\code{sample.phybreak} is deprecated).
#' 
#' @param x An object of class \code{phybreak}.
#' @param nsample The number of samples to take.
#' @param thin The thinning to use (values after every \code{thin}'th iteration will be included in the posterior). 
#'   Each iteration does one update of all parameters and tree updates with each host as focal host once.
#' @param classic The proportion of tree updates with the classic protocol (reference see below), in which within-host
#'   minitrees are proposed by simulating coalescent times and tree topology. In the current default protocol only 
#'   coalescent times are proposed with the minitree topology kept intact. This is followed by removing and reconnecting
#'   the sampling tips one by one. This results in better mixing of the mcmc-chain if there is much 
#'   genetic information (many SNPs) and/or if there are many possible within-host minitree topologies 
#'   (e.g. many samples per host). The classic protocol is faster in terms of updates/second and can thus be more efficient
#'   with little genetic information.
#' @param keepphylo The proportion of tree updates keeping the phylogenetic tree intact, only possible if there is one 
#'   sample per host and the \code{wh.model = "linear"} with complete bottleneck.
#' @param withinhost_only The proportion of tree updates in which only the within-host minitree is sampled, and 
#'   the transmission tree and infection times are kept unchanged.
#' @param parameter_frequency The relative frequency by which the model parameters are updated relative to updating each host.
#' @param status_interval The number of seconds between each on-screen print of the progress of the mcmc-chain.
#' @return The \code{phybreak} object used to call the function, including (additional) samples from the posterior.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak-object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 20)
#' MCMCstate <- sample_phybreak(MCMCstate, nsample = 50, thin = 2)
#' @export
sample_phybreak_parallel <- function(x, nsample, thin = 1, thinswap = 1, classic = 0, keepphylo = 0, withinhost_only = 0, 
                            parameter_frequency = 1, status_interval = 10,
                            histtime = -1e5, history = FALSE,
                            nchains = 1, heats = NULL, all_chains = FALSE, 
                            outfile = "") {
  
  if(!("parallel" %in% .packages(TRUE))) {
    stop("package 'parallel' should be installed for this function")
  }
  if(!("parallel" %in% .packages(FALSE))) {
    stop("package 'parallel' is not attached")
  }
  
  if(!("Rdsm" %in% .packages(TRUE))) {
    stop("package 'Rdsm' should be installed for this function")
  }
  if(!("Rdsm" %in% .packages(FALSE))) {
    stop("package 'Rdsm' is not attached")
  }
  
  ### tests
  if(nsample < 1) stop("nsample should be positive")
  if(thin < 1) stop("thin should be positive")
  if(is.null(x$p$wh.bottleneck)) {
    x$p$wh.bottleneck <- choose_whbottleneck("auto", x$p$wh.model)
  }
  if(classic < 0 | classic > 1) stop("classic should be a fraction")
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(withinhost_only < 0 | withinhost_only > 1) stop("withinhost_only should be a fraction")
  if(withinhost_only + keepphylo + classic > 1) stop("classic + keepphylo + withinhost_only should be a fraction")
  if(keepphylo > 0) {
    if(any(duplicated(x$d$hostnames)) || !(x$p$wh.model %in% c(3, "linear")) || x$p$wh.bottleneck == "wide") {
      keepphylo <- 0
      warning("model incompatible with keepphylo-updates: they will not be used", immediate. = TRUE)
    } 
  }
  if(nchains==1){
   return(sample_phybreak(x, nsample, thin = 1, thinswap = 1, classic = 0, keepphylo = 0, withinhost_only = 0, 
                          parameter_frequency = 1, status_interval = 10,
                          histtime = -1e5, history = FALSE))
  }
  
  ### add distance model if not present
  if(is.null(x$p$dist.model)) {
    x$p$dist.model <- "none"
    x$p$dist.exponent <- 2
    x$p$dist.scale <- 1
    x$p$dist.mean <- 1
    x$s$dist.e <- c()
    x$s$dist.s <- c()
    x$s$dist.m <- c()
    x$h$est.dist.s <- FALSE
    x$h$est.dist.e <- FALSE
    x$h$est.dist.m <- FALSE
  }
  
  protocoldistribution <- c(1 - classic - keepphylo - withinhost_only, classic, keepphylo, withinhost_only)
  
  ### create room in s to add the new posterior samples
  s.post <- list(inftimes = with(x, cbind(s$inftimes, matrix(NA, nrow = p$obs, ncol = nsample))),
                  infectors = with(x, cbind(s$infectors, matrix(NA, nrow = p$obs, ncol = nsample))),
                  nodetimes = with(x, cbind(s$nodetimes, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                  nodehosts = with(x, cbind(s$nodehosts, matrix(NA, nrow = d$nsamples - 1, ncol = nsample))), 
                  nodeparents = with(x, cbind(s$nodeparents, matrix(NA, nrow = 2 * d$nsamples - 1, ncol = nsample))),
                  introductions = c(sum(x$s$infectors==0), rep(NA, nsample - 1)),
                  mu = c(x$s$mu, rep(NA, nsample)), 
                  hist_dens = c(x$s$hist_dens, rep(NA, nsample)),
                  mG = c(x$s$mG, rep(NA, nsample)), 
                  mS = c(x$s$mS, rep(NA, nsample)), 
                  tG = c(x$s$tG, rep(NA, nsample)),
                  tS = c(x$s$tS, rep(NA, nsample)),
                  wh.h = c(x$s$wh.h, rep(NA, nsample)), 
                  wh.s = c(x$s$wh.s, rep(NA, nsample)), 
                  wh.e = c(x$s$wh.e, rep(NA, nsample)), 
                  wh.0 = c(x$s$wh.0, rep(NA, nsample)), 
                  dist.e = c(x$s$dist.e, rep(NA, nsample)), 
                  dist.s = c(x$s$dist.s, rep(NA, nsample)), 
                  dist.m = c(x$s$dist.m, rep(NA, nsample)), 
                  logLik = c(x$s$logLik, rep(NA, nsample)),
                  historyinf = c(x$s$historyinf, rep(NA, nsample)),
                  heat = c(x$s$heat, rep(NA, nsample)))
  
  ### Set heats for MC3, heats = 1 if there is one chain
  if (is.null(heats))
    heats <- 1/(1+1*(1:nchains-1))
  else if (inherits(heats, "numeric") & length(heats) != nchains)
    stop("length of heats is not the same as number of chains")
  else if (!inherits(heats, "numeric"))
    stop("heats is not a numeric vector")

  nswaps <- nsample / thinswap
  
  ### Set up cluster
  cores <- detectCores()
  cl <- makeCluster(min(nchains, cores-1), outfile = outfile)
  mgrinit(cl)
  makebarr(cl)
  
  ### Make shared variable
  mgrmakevar(cl, "shared_heats", 1, nchains)
  shared_heats[,] <- heats
  
  mgrmakevar(cl, "shared_lik", 1, nchains)
  
  ### Export variables and functions
  clusterExport(cl, varlist = c("x", "nsample", "nswaps", "thin", "thinswap", "protocoldistribution",
                                "s.post", "status_interval", "parameter_frequency",
                                "history"), 
                envir = environment())
  
  if (outfile != "")
    print(paste("Intermediate results of MCMC chain with heat beta=1 are shown in file:", outfile))
  
  ### Process MCMC for each chain in parallel and share heats and likelihoods
  process_mcmc <- function(shared_heats, shared_lik) {
    
    me <- myinfo$id
    if (nswaps == 1) { 
      if(me == 1)
        return(sample_phybreak(x, nsample = thinswap, heats=shared_heats[1,me],
                               history = history, keep_history = FALSE))
      else
        return(sample_phybreak(x, nsample = thinswap, heats=shared_heats[1,me],
                               history = history, keep_history = FALSE, 
                               verbose = 0))
    } else {
      if (me == 1)
        s <- sample_phybreak(x, nsample = thinswap, heats=shared_heats[1,me],
                             history = history, keep_history = TRUE)
      else
        s <- sample_phybreak(x, nsample = thinswap, heats=shared_heats[1,me],
                             history = history, keep_history = TRUE, 
                             verbose = 0)
    }
    if (nswaps == 2) {
      return(sample_phybreak(s, nsample = thinswap, heats=shared_heats[1,me],
                             history = TRUE, 
                             verbose = 0))
    } else {
      for (i in 2:(nswaps-1)){
        if (shared_heats[1,me] == 1)
          s <- sample_phybreak(s, nsample = thinswap, heats=shared_heats[1,me],
                               history = TRUE, keep_history = TRUE, 
                               verbose = 2)
        else 
          s <- sample_phybreak(s, nsample = thinswap, heats=shared_heats[1,me],
                               history = TRUE, keep_history = TRUE, 
                               verbose = 0)
        shared_lik[1, me] <- s$s$logLik[length(s$s$logLik)]
        #print(shared_lik[,])
        barr()
        if (me == 1){
          shared_heats[,] <- phybreak:::swap_heats(shared_heats[,], shared_lik[,])
        }
        barr()
      }
      
      s <- sample_phybreak(s, nsample = thinswap, heats=shared_heats[1,me],
                           history = TRUE, verbose = 0)
      
      print(warnings())
      return(s)
    }
  }
  
  clusterExport(cl, "process_mcmc", envir = environment())
  clusterEvalQ(cl, library(phybreak))
  clusterEvalQ(cl, library(Rdsm))

  posts <- clusterEvalQ(cl, process_mcmc(shared_heats, shared_lik))
  stopCluster(cl)

  ### Sort samples per chain
  heat_index <- matrix(nrow = nsample, ncol = nchains)
  for (i in 1:nchains) {
    for (j in 1:length(posts)) {
      heat_index[which(posts[[j]]$s$heat == heats[i]),j] <- i
    }
  }

  # e <- posts[[which(heat_index[nsample,] == 1)]]$p
  # for (i in ls(envir=e)) copy2pbe0(i, e)

  s.posts <- lapply(1:length(heats), function(nheat){
    s <- s.post
    for (chain in 1:nchains){
      smp <- which(posts[[chain]]$s$heat == heats[nheat])
      for (n in names(s)){
        if(inherits(s[[n]], "matrix"))
          s[[n]][,smp] <- posts[[chain]]$s[[n]][,smp]
        else
          s[[n]][smp] <- posts[[chain]]$s[[n]][smp]
      }
    }
    s$chain <- posts[[nheat]]$s$heat
    return(s)
  })

  #res <- s.posts[[1]]

  res <- posts[[which(heat_index[nsample,]==1)]]
  res$s <- s.posts[[1]]

  if(all_chains)
    return(list(phyb.obj = res, chains = s.posts))
  else
    return(res)
}
  
  
  