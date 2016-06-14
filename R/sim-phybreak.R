### simulation of outbreaks in obkData format, according to the phybreak-model ###

### calls:
# .sim.outbreak.size
# .sim.outbreak
# .sim.phylotree
# .sim.sequences
# .makephylo.phybreak
### required packages:
# OutbreakTools
# phangorn
sim.phybreak <- function(obsize = 50, popsize = NA, 
                         R0 = 1.5, shape.gen = 10, mean.gen = 1, 
                         shape.sample = 10, mean.sample = 1,
                         wh.model = 3, slope = 1, 
                         mutrate = 0.0001, sequence.length = 10000) {
  ### tests
  if(all(is.na(c(obsize,popsize)))) {
    stop("give an outbreak size (obsize) and/or a population size (popsize)")
  }
  if(all(!is.na(c(obsize,popsize))) && obsize > popsize) {
    stop("outbreak size (obsize) cannot be larger than population size (popsize)")
  }
  
  ### simulate step by step
  if(is.na(obsize)) {
     res <- .sim.outbreak(popsize, R0, shape.gen, mean.gen,
                                      shape.sample, mean.sample)
     obsize <- res$obs
     if(obsize == 1) return(c("Outbreak size = 1"))
  } else {
    res <- .sim.outbreak.size(obsize, popsize, R0, shape.gen, mean.gen,
                              shape.sample, mean.sample)
  }
  res <- .sim.phylotree(res, wh.model, slope)
  res <- .sim.sequences(res, mutrate, sequence.length)

  ### make an obkData object
  treesout <- vector('list',1)
  treesout[[1]] <- .makephylo.phybreak(res$nodetimes, res$nodeparents, 1:obsize)
  class(treesout) <- "multiPhylo"
  
  toreturn <- new("obkData",
                  individuals = data.frame(
                    infector = res$infectors,
                    date = as.Date(res$infectiontimes, origin = "2000-01-01"),
                    row.names = 1:res$obs),
                  dna = list(SNPs = res$SNPlist), dna.date = as.Date(res$sampletimes, origin = "2000-01-01"),
                  dna.individualID = 1:res$obs, trees = treesout)
  return(toreturn)
}
