### phybreakdata class constructor ###

### phybreak functions called ###
# .rinftimes(sampling times, sample.mean, sample.shape)
# .rinfectors(infection times, generation.mean, generation.shape)
# .samplecoaltimes(tip times, wh.model, wh.slope)
# .sampletopology(nodes, nodetimes, nodetypes, rootnode, wh.model)
# .distmatrix(SNP, SNPfrequencies)

#' Create a phybreakdata-object from raw data.
#' 
#' phybreakdata takes as data sequences and sampling times and makes a \code{phybreakdata} object.
#' Each sample is assumed to be associated with a separate host. The number of sequences
#' should be equal to the length of the sampling time vector, and the position identifies the host (unless named
#' vectors are provided). Host names can be provided separately;
#' otherwise it will be tried to extract them from the sequences or sampling times.
#' It is also possible to include (otherwise unobserved) simulated data: sim.infection times, sim.infectors, 
#' and a (phylogenetic) sim.tree. This is done automatically when using \code{\link{sim.phybreak}}. 
#' 
#' @param sequences Sequence data of class \code{'DNAbin'}, \code{'phyDat'}, or a \code{matrix} with nucleotides, 
#'  each row a host, each column a nucleotide). All nucleotides that are not \code{'a'}, \code{'c'}, \code{'g'}, 
#'  or \code{'t'}, will be turned into \code{'n'}. 
#' @param sample.times A vector of sampling times (\code{numerical} or \code{Date}).
#' @param host.names A vector with host names.
#' @param sim.infection.times A vector with infection times (\code{numerical} or \code{Date}).
#' @param sim.infectors A vector with infectors, either by name or by position (use 0 for the index case).
#' @param sim.tree A tree of class \code{'phylo'}, with tip names identifying the hosts.
#' @return An object of class \code{phybreakdata} with the following elements
#'   \describe{
#'     \item{sequences}{a \code{'phyDat'}-object with the sequence data.}
#'     \item{sample.times}{a named \code{vector} with the sample times.}
#'     \item{sim.infection.times}{a named \code{vector} with the (simulated) infection times (if provided).}
#'     \item{sim.infectors}{a named \code{vector} with the (simulated) infectors (if provided).}
#'     \item{sim.tree}{a \code{'phylo'}-object with the (simulated) phylogenetic tree (if provided).}
#'   }
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1101/069195}{Klinkenberg et al, on biorXiv}.
#' @examples 
#' sampletimedata <- c(0,2,2,4,4)
#' sampleSNPdata <- matrix(c("a","a","a","a","a",
#'                           "a","c","c","c","c",
#'                           "t","t","t","g","g"), nrow = 5, 
#'                           dimnames = list(LETTERS[1:5], NULL))
#' dataset <- phybreakdata(sequences = sampleSNPdata, sample.times = sampletimedata)
#' @export
phybreakdata <- function(sequences, sample.times, host.names = NULL, sim.infection.times = NULL,
                         sim.infectors = NULL, sim.tree = NULL) {
  
  ##########################################################
  ### testing the input: first the essential information ###
  ##########################################################
  if(!any(class(sequences) %in% c("DNAbin", "phyDat", "matrix"))) {
    stop("sequences should be of class \"DNAbin\", \"phyDat\", or \"character\"")
  }
  if(inherits(sequences, "matrix") && !inherits(sequences[1], "character")) {
    stop("sequences matrix should contain \"character\" elements")
  }
  if(inherits(sequences, c("DNAbin", "phyDat"))) {
    sequences <- as.character(sequences)
  }
  if(!inherits(sample.times, c("Date", "numeric", "integer"))) {
    stop("sample.times should be numeric or of class \"Date\"")
  }
  if(nrow(sequences) != length(sample.times)) {
    stop("numbers of sequences and sample.times don't match")
  }
  if(!is.null(host.names)) {
    if(length(host.names) != length(sample.times)) {
      stop("length of host.names does not match number of sequences")
    }
    if(is.null(names(sample.times))) {
      names(sample.times) <- host.names
    } else if (all(names(sample.times) %in% host.names)) {
      sample.times <- sample.times[host.names]
    } else {
      warning("names in sample.times don't match host.names and are therefore overwritten")
      names(sample.times) <- host.names
    }
    if(is.null(row.names(sequences))) {
      row.names(sequences) <- host.names
    } else if (all(row.names(sequences) %in% host.names)) {
      sequences <- sequences[host.names, ]
    } else {
      warning("names in sequences don't match host.names and are therefore overwritten")
      row.names(sequences) <- host.names
    }
  } else if(!is.null(row.names(sequences))) {
    host.names <- row.names(sequences)
    if(is.null(names(sample.times))) {
      names(sample.times) <- host.names
    } else if (all(names(sample.times) %in% host.names)) {
      sample.times <- sample.times[host.names]
    } else {
      warning("names in sample.times don't match sequence names and are therefore overwritten")
      names(sample.times) <- host.names
    }
  } else if(!is.null(names(sample.times))) {
    host.names <- names(sample.times)
    row.names(sequences) <- host.names
  } else {
    host.names <- paste0("host.", 1:length(sample.times))
    row.names(sequences) <- host.names
    names(sample.times) <- host.names
  }
  
  #################################################
  ### place essential information in outputlist ###
  #################################################
  if(length(setdiff(sequences,c("a","c","g","t")))) warning("all nucleotides other than actg are turned into n")
  sequences[sequences != "a" & sequences != "c" & sequences != "g" & sequences != "t"] <- "n"
  sequences <- phangorn::as.phyDat(sequences)
  
  res <- list(
    sequences = sequences,
    sample.times = sample.times
  )
  ######################################################################
  ### testing and adding the input: then the rest of the information ###
  ######################################################################
  if(!is.null(sim.infection.times)) {
    if(class(sim.infection.times) != class(sample.times)) {
      stop("sim.infection.times should be of same class as sample.times")
    }
    if(is.null(names(sim.infection.times))) {
      names(sim.infection.times) <- host.names
    } else if (all(names(sim.infection.times) %in% host.names)) {
      sim.infection.times <- sim.infection.times[host.names]
    } else {
      warning("names in sim.infection.times don't match host.names and are therefore overwritten")
      names(sim.infection.times) <- host.names
    }
    if(!all(sim.infection.times < sample.times)) {
      stop("all infection times should be before the sampling times")
    }
    if(length(host.names) != length(sim.infection.times)) {
      stop("length of sim.infection.times does not match number of hosts")
    }
    res <- c(res, list(sim.infection.times = sim.infection.times))
  }
  if(!is.null(sim.infectors)) {
    if(!inherits(sim.infectors, c("character", "numeric", "integer"))) {
      stop("sim.infectors should be numeric (referring to position of infector) or character (referring to host names)")
    }
    if(inherits(sim.infectors, c("numeric", "integer")) && !all(sim.infectors %in% 0:length(host.names))) {
      stop("sim.infectors should be integers between 0 and number of hosts")
    }
    if(inherits(sim.infectors, c("character")) && !all(sim.infectors %in% c(0, "index", host.names))) {
      stop("not all sim.infectors are proper host names")
    }
    if(inherits(sim.infectors, c("character")) && any(sim.infectors == "0")) {
      sim.infectors[sim.infectors == "0"] <- "index"
    }
    if(length(sim.infectors) != length(sample.times)) {
      stop("length of sim.infectors does not match number of hosts")
    }
    if(is.null(names(sim.infectors))) {
      names(sim.infectors) <- host.names
    } else if (all(names(sim.infectors) %in% host.names)) {
      sim.infectors <- sim.infectors[host.names]
    } else {
      warning("names in sim.infectors don't match host.names and are therefore overwritten")
      names(sim.infectors) <- host.names
    }
    res <- c(res, list(sim.infectors = sim.infectors))
  }
  if(!is.null(sim.tree)) {
    if(!inherits(sim.tree, "phylo")) {
      stop("sim.tree should be of class phylo")
    }
    if(!all(sim.tree$tip.label %in% host.names)) {
      stop("names in sim.tree don't match host names")
    }
    res <- c(res, list(sim.tree = sim.tree))
  }
  
  class(res) <- "phybreakdata"
  return(res)
}



