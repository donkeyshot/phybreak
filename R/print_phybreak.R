#' Description of phybreak object.
#' 
#' Give a short description of the data and the model, with the current parameter values.
#' 
#' @param x A \code{phybreak} object.
#' @param printlen The number of host names to print.
#' @param ... Some methods for this generic require additional arguments. None are used in this method.
#' @return NULL
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' simulation <- sim_phybreak()
#' x <- phybreak(simulation)
#' print(x)
#' @export
print.phybreak <- function(x, ...) {
  print_phybreak_general(x)
  cat("\n\n")
  print_data(x)
  cat("\n\n")
  print_model(x)
}

print_phybreak_general <- function(x) {
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  cat("A phybreak object with ", length(x$s$mu), " posterior samples.", sep = "")
}

print_model <- function(x) {
  if (!inherits(x, "phybreak")) {
    stop("object must be of class \"phybreak\"")
  }
  
  cat("The model:\n")
  cat("* generation interval ~ Gamma(shape = ", signif(x$p$gen.shape, 3), ", mean = ", signif(x$p$gen.mean),")\n", sep = "")
  if(x$h$est.mG) {
    cat("\tmean estimated with prior mean = ", x$h$mG.av, " and sd = ", x$h$mG.sd, "\n", sep = "")
  } else {
    cat("\tmean not estimated\n")
  }
  cat("* sampling interval ~ Gamma(shape = ", signif(x$p$sample.shape, 3), ", mean = ", signif(x$p$sample.mean),")\n", sep = "")
  if(x$h$est.mS) {
    cat("\tmean estimated with prior mean = ", x$h$mS.av, " and sd = ", x$h$mS.sd, "\n", sep = "")
  } else {
    cat("\tmean not estimated\n")
  }
  cat("* within-host model is ")
  switch(x$p$wh.model,
         single = "\"single\"\n",
         infinite = "\"infinite\"\n",
         linear = {
           cat("\"linear\" with slope = ", signif(x$p$wh.slope, 3), "\n\tslope ", sep = "")
           if(x$h$est.wh.s) cat("estimated with Gamma prior (shape = ", x$h$wh.s.sh,
                                ", mean = ", x$h$wh.s.av, ")\n", sep = "") else cat("not estimated\n") 
         },
         exponential = {
           cat("\"exponential\" with initial level = ", signif(x$p$wh.level, 3), " and exponent = ", 
               signif(x$p$wh.exponent, 3), "\n\tlevel ", sep = "")
           if(x$h$est.wh.0) cat("estimated with Gamma prior (shape = ", x$h$wh.0.sh,
                                ", mean = ", x$h$wh.0.av, ")\n", sep = "") else cat("not estimated\n") 
           cat("\texponent ")
           if(x$h$est.wh.e) cat("estimated with Gamma prior (shape = ", x$h$wh.e.sh,
                                ", mean = ", x$h$wh.e.av, ")\n", sep = "") else cat("not estimated\n") 
         },
         constant = {
           cat("\"constant\" with level = ", signif(x$p$wh.level, 3), "\n\tlevel ", sep = "")
           if(x$h$est.wh.0) cat("estimated with Gamma prior (shape = ", x$h$wh.0.sh,
                                ", mean = ", x$h$wh.0.av, ")\n", sep = "") else cat("not estimated\n") 
         })
  cat("* mutation rate mu = ", signif(x$p$mu, 3), ", estimated with log-uniform prior", sep = "")
  
}

print_data <- function(x, printlen = 6, ...) {
  cat("The data:\n")
  cat("* an outbreak of ", x$p$obs, " hosts, named ", sep = "")
  if(printlen >= x$p$obs) {
    cat(paste(x$d$hostnames[1:x$p$obs], collapse = ", "), "\n")
  } else {
    cat(paste(x$d$hostnames[1:printlen], collapse = ", "), ", ...", "\n", sep = "")
  }
  if(inherits(x$d$sample.times, "Date")) {
    cat("* the outbreak started ", as.character(min(x$d$sample.times)), 
        " and ended ", as.character(max(x$d$sample.times)), "\n", sep = "")
  } else {
    cat("* the last sample was taken ", max(x$d$sample.times) - min(x$d$sample.times), " time units after the first\n", sep = "")
  }
  cat("* the data contain ", x$d$nsamples, " sequences, with a total of ", x$d$nSNPs, " SNPs", sep = "")
}

