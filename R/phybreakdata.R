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
#' If no host names are provided, each sample is assumed to be associated with a separate host. The number of sequences
#' should be equal to the length of the sampling time vector, and the position identifies the host (unless named
#' vectors are provided). Sample names can be provided separately;
#' otherwise it will be tried to extract them from the sequences or sampling times.
#' It is also possible to include (otherwise unobserved) simulated data: sim.infection times, sim.infectors, 
#' and a (phylogenetic) sim.tree. This is done automatically when using \code{\link{sim.phybreak}}. 
#' 
#' @param sequences Sequence data of class \code{'DNAbin'}, \code{'phyDat'}, or a \code{matrix} with nucleotides, 
#'  each row a host, each column a nucleotide). In a matrix, nucleotides should be lower-case letters. All undefined 
#'  nucleotides or ambiguity codes will be turned into \code{'n'}. 
#' @param sample.times A vector of sampling times (\code{numerical} or \code{Date}).
#' @param sample.names A vector with sample names.
#' @param host.names A vector with host names. The vector identifies the host for each sample, so should be of the same
#'  length as \code{sample.times}. 
#' @param sim.infection.times A vector with infection times (\code{numerical} or \code{Date}).
#' @param sim.infectors A vector with infectors, either by name or by position (use 0 for the index case).
#' @param sim.tree A tree of class \code{'phylo'}, with tip names identifying the hosts.
#' @return An object of class \code{phybreakdata} with the following elements
#'   \describe{
#'     \item{sequences}{a \code{'phyDat'}-object with the sequence data.}
#'     \item{sample.times}{a named \code{vector} with the sample times.}
#'     \item{sample.hosts}{a named \code{vector} with the hosts from whom the samples have been taken.}
#'     \item{sim.infection.times}{a named \code{vector} with the (simulated) infection times (if provided).}
#'     \item{sim.infectors}{a named \code{vector} with the (simulated) infectors (if provided).}
#'     \item{sim.tree}{a \code{'phylo'}-object with the (simulated) phylogenetic tree (if provided).}
#'   }
#' @author Don Klinkenberg \email{don@@xs4all.nl}
#' @references \href{http://dx.doi.org/10.1371/journal.pcbi.1005495}{Klinkenberg et al. (2017)} Simultaneous 
#'   inference of phylogenetic and transmission trees in infectious disease outbreaks. 
#'   \emph{PLoS Comput Biol}, \strong{13}(5): e1005495.
#' @examples 
#' sampletimedata <- c(0,2,2,4,4)
#' sampleSNPdata <- matrix(c("a","a","a","a","a",
#'                           "a","c","c","c","c",
#'                           "t","t","t","g","g"), nrow = 5, 
#'                           dimnames = list(LETTERS[1:5], NULL))
#' dataset <- phybreakdata(sequences = sampleSNPdata, sample.times = sampletimedata)
#' @export
phybreakdata <- function(sequences, sample.times, sample.names = NULL, host.names = sample.names, 
                         sim.infection.times = NULL, sim.infectors = NULL, sim.tree = NULL) {
  
  ##########################################################
  ### testing the input: first the essential information ###
  ##########################################################
  if(!any(class(sequences) %in% c("DNAbin", "phyDat", "matrix"))) {
    stop("sequences should be of class \"DNAbin\", \"phyDat\", or \"character\"")
  }
  if(inherits(sequences, "matrix") && !inherits(sequences[1], "character")) {
    stop("sequences matrix should contain \"character\" elements")
  }
  if(inherits(sequences, c("DNAbin"))) {
    sequences <- phangorn::as.phyDat(sequences)
  }
  if(inherits(sequences, c("phyDat"))) {
    sequences <- as.character(sequences)
  }
  if(!inherits(sample.times, c("Date", "numeric", "integer"))) {
    stop("sample.times should be numeric or of class \"Date\"")
  }
  if(nrow(sequences) != length(sample.times)) {
    stop("numbers of sequences and sample.times don't match")
  }
  if(is.null(sample.names)) {
    if(!is.null(row.names(sequences))) {
      sample.names <- row.names(sequences)
    } else if(!is.null(names(sample.times))) {
      sample.names <- names(sample.times)
    } else if(!is.null(host.names) && length(sample.times) == length(unique(host.names))) {
    } else {
      sample.names <- paste0("sample.", 1:length(sample.times))
    }
  }
  if(length(sample.names) != length(sample.times)) {
    stop("length of sample.names does not match number of sequences")
  }
  if(is.null(row.names(sequences))) {
    row.names(sequences) <- sample.names
  } else if (all(row.names(sequences) %in% sample.names)) {
    sequences <- sequences[sample.names, ]
  } else {
    warning("names in sequences don't match sample.names and are therefore overwritten")
    row.names(sequences) <- sample.names
  }
  if(is.null(names(sample.times))) {
    names(sample.times) <- sample.names
  } else if (all(names(sample.times) %in% sample.names)) {
    sample.times <- sample.times[sample.names]
  } else {
    warning("names in sample.times don't match sample.names and are therefore overwritten")
    names(sample.times) <- sample.names
  }
  if(is.null(host.names)) host.names <- sample.names
  if(length(host.names) != length(sample.names)) {
    stop("length of host.names does not match number of sequences")
  }
  if(is.null(names(host.names))) {
    names(host.names) <- sample.names
  } else if(all(names(host.names) %in% sample.names)) {
    host.names <- host.names[sample.names]
  } else {
    warning("names in host.names don't match sample.names and are therefore overwritten")
    names(host.names) <- sample.names
  }

  ##########################################################################
  ### order the input host by host, hosts ordered by first sampling time ###
  ##########################################################################
  allhosts <- unique(host.names)
  allfirsttimes <- rep(FALSE, length(sample.times))
  sapply(allhosts, function(x) allfirsttimes[which(min(sample.times[host.names == x]) == sample.times)[1]] <<- TRUE)
  outputorderhosts <- order(sample.times[allfirsttimes])
  orderedhosts <- host.names[allfirsttimes][outputorderhosts]
  outputordersamples <- order(!allfirsttimes, match(host.names, orderedhosts), sample.times)
  sequences <- sequences[outputordersamples, ]
  sample.times <- sample.times[outputordersamples]
  sample.names <- sample.names[outputordersamples]
  host.names <- host.names[outputordersamples]

  #################################################
  ### place essential information in outputlist ###
  #################################################
  # if(length(setdiff(sequences,c("a","c","g","t","u","m","r","w",
  #                               "s","y","k","v","h","d","b","n",
  #                               "?","-")))) warning("all undefined nucleotide codes are turned into n")
  # sequences[!(sequences %in% c("a","c","g","t","u","m","r","w",
  #                              "s","y","k","v","h","d","b","n",
  #                              "?","-"))] <- "n"
  if(length(setdiff(sequences,c("a","c","g","t","n")))) warning("all nucleotide ambiguity codes are turned into n")
  sequences[!(sequences %in% c("a","c","g","t"))] <- "n"
  sequences <- phangorn::as.phyDat(sequences)
  
  res <- list(
    sequences = sequences,
    sample.times = sample.times,
    sample.hosts = host.names
  )
  ######################################################################
  ### testing and adding the input: then the rest of the information ###
  ######################################################################
  if(!is.null(sim.infection.times)) {
    if(class(sim.infection.times) != class(sample.times)) {
      stop("sim.infection.times should be of same class as sample.times")
    }
    if(length(orderedhosts) != length(sim.infection.times)) {
      stop("length of sim.infection.times does not match number of hosts")
    }
    if(is.null(names(sim.infection.times))) {
      sim.infection.times <- sim.infection.times[outputorderhosts]
      names(sim.infection.times) <- orderedhosts
    } else if (all(names(sim.infection.times) %in% orderedhosts)) {
      sim.infection.times <- sim.infection.times[orderedhosts]
    } else {
      warning("names in sim.infection.times don't match host.names and are therefore overwritten")
      sim.infection.times <- sim.infection.times[outputorderhosts]
      names(sim.infection.times) <- orderedhosts
    }
    if(!all(sim.infection.times < sample.times[1:length(allhosts)])) {
      stop("all infection times should be before the sampling times")
    }
    res <- c(res, list(sim.infection.times = sim.infection.times))
  }
  if(!is.null(sim.infectors)) {
    if(!inherits(sim.infectors, c("character", "numeric", "integer"))) {
      stop("sim.infectors should be numeric (referring to position of infector) or character (referring to host names)")
    }
    if(inherits(sim.infectors, c("numeric", "integer")) && !all(sim.infectors %in% 0:length(orderedhosts))) {
      stop("sim.infectors should be integers between 0 and number of hosts")
    }
    if(inherits(sim.infectors, c("character")) && !all(sim.infectors %in% c(0, "index", orderedhosts))) {
      stop("not all sim.infectors are proper host names")
    }
    if(inherits(sim.infectors, c("character")) && any(sim.infectors == "0")) {
      sim.infectors[sim.infectors == "0"] <- "index"
    }
    if(length(sim.infectors) != length(orderedhosts)) {
      stop("length of sim.infectors does not match number of hosts")
    }
    if(is.null(names(sim.infectors))) {
      sim.infectors <- sim.infectors[outputorderhosts]
      names(sim.infectors) <- orderedhosts
    } else if (all(names(sim.infectors) %in% orderedhosts)) {
      sim.infectors <- sim.infectors[orderedhosts]
    } else {
      warning("names in sim.infectors don't match host.names and are therefore overwritten")
      sim.infectors <- sim.infectors[outputorderhosts]
      names(sim.infectors) <- orderedhosts
    }
    if(inherits(sim.infectors, c("numeric", "integer"))) {
      sim.infectors <- c("index", orderedhosts)[sim.infectors + 1]
    }
    names(sim.infectors) <- orderedhosts
    res <- c(res, list(sim.infectors = sim.infectors))
  }
  if(!is.null(sim.tree)) {
    if(!inherits(sim.tree, "phylo")) {
      stop("sim.tree should be of class phylo")
    }
    if(!all(sim.tree$tip.label %in% sample.names)) {
      stop("names in sim.tree don't match sample names")
    }
    currenttipsinedge <- match(1:length(sample.times), sim.tree$edge[,2])
    tipreorder <- match(sim.tree$tip.label, sample.names)
    sim.tree$edge[currenttipsinedge, 2] <- tipreorder
    sim.tree$tip.label <- sample.names
    res <- c(res, list(sim.tree = sim.tree))
  }
  
  class(res) <- "phybreakdata"
  return(res)
}



