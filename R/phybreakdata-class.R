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
#' Each sample is assumed to be associated with a separate host. Additional data can also be given, 
#' such as infection times, infectors, and a phylogenetic tree. Host names can be provided separately;
#' otherwise it will be tried to extract them from the sequences or sampling times. The number of sequences
#' should be equal to the length of all data vectors, and the position identifies the host (unless named
#' vectors are provided).
#' 
#' @param sequences Sequence data of class \code{'DNAbin'}, \code{'phyDat'}, or a \code{matrix} with nucleotides, 
#'  each row a host, each column a nucleotide). All nucleotides that are not \code{'a'}, \code{'c'}, \code{'g'}, 
#'  or \code{'t'}, will be turned into \code{'n'}. 
#' @param sample.times A vector of sampling times (\code{numerical} or \code{Date}).
#' @param host.names A vector with host names.
#' @param infection.times A vector with infection times (\code{numerical} or \code{Date}).
#' @param infectors A vector with infectors, either by name or by position (use 0 for the index case).
#' @param tree A tree of class \code{\link{ape::phylo}}, with tip names identifying the hosts.
#' @return An object of class \code{phybreakdata} with the following elements
#'   \describe{
#'     \item{sequences}{a \code{'phyDat'}-object with the sequence data.}
#'     \item{sample.times}{a named \code{vector} with the sample times.}
#'     \item{infection.times}{a named \code{vector} with the infection times (if provided).}
#'     \item{infectors}{a named \code{vector} with the infectors (if provided).}
#'     \item{tree}{a \code{'phylo'}-object with the phylogenetic tree (if provided).}
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
phybreakdata <- function(sequences, sample.times, host.names = NULL, infection.times = NULL,
                         infectors = NULL, tree = NULL) {
  
  ##########################################################
  ### testing the input: first the essential information ###
  ##########################################################
  if(!any(class(sequences) %in% c("DNAbin", "phyDat", "matrix"))) {
    stop("sequences should be of class \"DNAbin\", \"phyDat\", or \"character\"")
  }
  if(inherits(sequences, "matrix") && !inherits(sequences[1], "character")) {
    stop("sequences matrix should contain \"character\" elements")
  }
  if(inherits(sequences, c("DNAbin"), "phyDat")) {
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
  if(!is.null(infection.times)) {
    if(!inherits(infection.times, c("Date", "numeric", "integer"))) {
      stop("sample.times should be numeric or of class \"Date\"")
    }
    if(!all(infection.times < sample.times)) {
      stop("all infection times should be before the sampling times")
    }
    if(length(host.names) != length(infection.times)) {
      stop("length of infection.times does not match number of hosts")
    }
    if(is.null(names(infection.times))) {
      names(infection.times) <- host.names
    } else if (all(names(infection.times) %in% host.names)) {
      infection.times <- infection.times[host.names]
    } else {
      warning("names in infection.times don't match host.names and are therefore overwritten")
      names(infection.times) <- host.names
    }
    res <- c(res, list(infection.times = infection.times))
  }
  if(!is.null(infectors)) {
    if(!inherits(infectors, c("character", "numeric", "integer"))) {
      stop("infectors should be numeric (referring to position of infector) or character (referring to host names)")
    }
    if(inherits(infectors, c("numeric", "integer")) && !all(infectors %in% 0:length(host.names))) {
      stop("infectors should be integers between 0 and number of hosts")
    }
    if(inherits(infectors, c("character")) && !all(infectors %in% c(0, "index", host.names))) {
      stop("not all infectors are proper host names")
    }
    if(inherits(infectors, c("character")) && any(infectors == "0")) {
      infectors[infectors == "0"] <- "index"
    }
    if(length(infectors) != length(infection.times)) {
      stop("length of infectors does not match number of hosts")
    }
    if(is.null(names(infectors))) {
      names(infectors) <- host.names
    } else if (all(names(infectors) %in% host.names)) {
      infectors <- infectors[host.names]
    } else {
      warning("names in infectors don't match host.names and are therefore overwritten")
      names(infectors) <- host.names
    }
    res <- c(res, list(infectors = infectors))
  }
  if(!is.null(tree)) {
    if(!inherits(tree, "phylo")) {
      stop("tree should be of class phylo")
    }
    if(!all(tree$tip.label %in% host.names)) {
      stop("names in tree don't match host names")
    }
    res <- c(res, list(tree = tree))
  }
  
  class(res) <- "phybreakdata"
  return(res)
}



