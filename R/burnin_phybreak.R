#' MCMC updating of a phybreak-object.
#' 
#' This function allows the MCMC chain to burn in. If used after samples have been taken (with \code{\link{sample_phybreak}}), 
#'   these samples will be returned unchanged in the output (\code{burnin.phybreak} is deprecated).
#' 
#' @param x An object of class \code{phybreak}.
#' @param ncycles Number of iterations to be carried out. Each iteration does one update of all parameters and
#'   tree updates with each host as focal host once.
#' @param keepphylo The proportion of tree updates keeping the phylotree intact, only possible if there is one sample per host
#'   and the \code{wh.model = "linear"}. If \code{NULL} (default), it is set to 0.2 only in that case, otherwise to 0.
#' @param withinhost_only The proportion of tree updates in which only the within-host minitree is sampled, and 
#'   the transmission tree and infection times are kept unchanged.
#' @return The \code{phybreak} object provided as input, with variables and parameters changed due to the updating.
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' #First create a phybreak object
#' simulation <- sim_phybreak(obsize = 5)
#' MCMCstate <- phybreak(dataset = simulation)
#' 
#' MCMCstate <- burnin_phybreak(MCMCstate, ncycles = 50)
#' @export
burnin_phybreak <- function(x, ncycles, keepphylo = NULL, withinhost_only = 0) {
  ### tests
  if(ncycles < 1) stop("ncycles should be positive")
  if(is.null(keepphylo)) {
    keepphylo <- 0.2
  }
  if(keepphylo > 0) {
    if(any(duplicated(x$d$hostnames)) || !(x$p$wh.model %in% c(3, "linear")) ) {
      keepphylo <- 0
      message("keepphylo = 0")
    } else {
      message(paste0("keepphylo = ", keepphylo))
    }
  }
  if(keepphylo < 0 | keepphylo > 1) stop("keepphylo should be a fraction")
  if(withinhost_only < 0 | withinhost_only > 1) stop("withinhost_only should be a fraction")
  if(withinhost_only + keepphylo > 1) stop("keepphylo + withinhost_only should be a fraction")
  
  build_pbe(x)
  
  message(paste0("   cycle      logLik         mu  gen.mean  sam.mean parsimony (nSNPs = ", pbe0$d$nSNPs, ")"))
  print_screen_log(0)
  
  curtime <- Sys.time()
tostop <<- FALSE
  for (rep in 1:ncycles) {
    if(Sys.time() - curtime > 10) {
      print_screen_log(rep)
      curtime <- Sys.time()
    }
# cat(paste0("rep ", rep, ": "))
    for (i in sample(x$p$obs)) {
      # cat(paste0(i, " "))
      # if(rep == 3 && i == 8) {Sys.sleep(1); tostop <<- TRUE}
      if (runif(1) < 1 - keepphylo - withinhost_only) 
        update_host(i) else  if (runif(1) < keepphylo/(keepphylo + withinhost_only)) {
          update_host_keepphylo(i)
        } else update_host_withinhost(i)
    }
    if (x$h$est.mG)
      update_mG()
    if (x$h$est.mS)
      update_mS()
    if (x$h$est.wh.s)
      update_wh_slope()
    if (x$h$est.wh.e)
      update_wh_exponent()
    if (x$h$est.wh.0)
      update_wh_level()
    update_mu()
  }
  
  res <- destroy_pbe(x$s)
  
  return(res)
}


#' @rdname burnin_phybreak
#' @export
burnin.phybreak <- function(...) {
  .Deprecated("burnin_phybreak")
  burnin_phybreak(...)
}

print_screen_log <- function(iteration) {
  message(paste0(
    stringr::str_pad(iteration, 8),
    stringr::str_pad(round(pbe0$logLikseq + pbe0$logLiksam + pbe0$logLikgen + pbe0$logLikcoal, 2), 12),
    stringr::str_pad(signif(pbe0$p$mu, 3), 11),
    stringr::str_pad(signif(pbe0$p$gen.mean, 3), 10),
    stringr::str_pad(signif(pbe0$p$sample.mean, 3), 10),
    stringr::str_pad(phangorn::parsimony(
      phybreak2phylo(environment2phybreak(pbe0$v)), pbe0$d$sequences), 10)))
}
