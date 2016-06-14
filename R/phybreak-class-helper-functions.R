### functions used by phybreak class constructor ###

### random infection times given sampling times and sampling interval distribution
### called from: 
# phybreak
.rinftimes <- function(st, meanS, shapeS) {
  ### tests
  if(class(st) != "numeric" && class(st) != "integer") {
    stop(".rinftimes called with non-numeric sampling times")
  }
  if(meanS <= 0) stop(".rinftimes called with non-positive mean sampling interval")
  if(shapeS <= 0) stop(".rinftimes called with non-positive shape parameter")
  
  ### function body
  st - rgamma(length(st), shape = shapeS, scale = meanS/shapeS)
}

### random infectors times given infection times and generation interval distribution
### called from:
# phybreak
.rinfectors <- function(it, meanG, shapeG) {
  ### tests
  if(class(it) != "numeric" && class(it) != "integer") {
    stop(".rinfectors called with non-numeric infection times")
  }
  if(sum(it == min(it)) > 1) stop("rinfectors with >1 index case")
  if(meanG <= 0) stop(".rinfectors called with non-positive mean generation interval")
  if(shapeG <= 0) stop(".rinfectors called with non-positive shape parameter")
  
  ### function body
  res <- rep(0,length(it))
  for(i in 1:length(it)) {
    if(it[i] > min(it)) {
      dist <- dgamma(it[i] - it, shape = shapeG, scale = meanG/shapeG)
      dist[i] <- 0
      res[i] <- sample(length(it), 1, prob = dist)
    }
  }
  return(res)
}

### distance matrix between sequences given SNP data
### called from:
# phybreak
.distmatrix <- function(SNPs, SNPfreqs) {
  ### tests
  if(ncol(SNPs) != length(SNPfreqs)) {
    stop(".distmatrix called with different SNP numbers in SNPs and SNPfreqs")
  }
  
  ### function body
  res <- matrix(0, nrow = nrow(SNPs), ncol = nrow(SNPs))
  #count SNPs excluding "n"
  for(i in 1:nrow(SNPs)) {
    for(j in i:nrow(SNPs)) {
      res[i,j] <- sum((SNPs[i,]!=SNPs[j,] & SNPs[i,]!="n" & SNPs[j,]!="n")*SNPfreqs)
      res[j,i] <- res[i,j]
    }
  }
  #prob of SNP per nucleotide in most distant entry
  nscore <- max(res)/sum(SNPfreqs)
  #add nscore for each missing nucleotide
  for(i in 1:nrow(SNPs)) {
    for(j in i:nrow(SNPs)) {
      res[i,j] <- res[i,j] + sum((SNPs[i,]=="n" | SNPs[j,]=="n")*SNPfreqs)*nscore
      res[j,i] <- res[i,j]
    }
  }
  
  #add 1 to avoid division by 0, and make distances proportional
  return((res+1)/max(res+1))
}
