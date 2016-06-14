### function to remove the transmission nodes from a nodeparents-set of a phybreak-object
### called by:
# .makephyloparsets
.makephyloparset <- function(parentset) {
  obs <- (1 + length(parentset))/3
  res <- parentset
  while(max(res) >= 2*obs) {
    res[res >= 2*obs] <- parentset[res][res >= 2*obs]
  }
  return(res[1 : (2*obs-1)])
}

### function to remove the transmission nodes from a matrix with nodeparents-sets
### called by:
# .phylotree
### calls:
# .makephyloparset
.makephyloparsets <- function(nodeparentsets) {
  apply(nodeparentsets,
        MARGIN = 2, .makephyloparset)
}


