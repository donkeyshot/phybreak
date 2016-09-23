#include <Rcpp.h>
#include <vector>
#include <string>
using namespace Rcpp;

// Calculation of log-likelihood of sequence data
// given a phylo-transmission tree

// [[Rcpp::export(name=".likseq")]]
double likseq(IntegerVector SNPs, IntegerVector SNPfreqs,
              IntegerVector nodeparents, NumericVector nodetimes,
              double mutrate, int obs) {
  int nnodes = 2*obs - 1;
  int nSNPs = SNPs.size()/obs;
  NumericVector likarray(nnodes * nSNPs * 4);
  NumericVector nlens(nodetimes.size());
  NumericVector SNPsums(nSNPs);
  LogicalVector routefree(nodetimes.size());
  int curnode;
  int nextnode;
  int rootnode = 0;
  double edgelen;
  double totprob;
  double result;
  
  for(int i = 0; i < obs; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      if(SNPs[i*nSNPs + j] == 1) {
        likarray[(i * nSNPs + j) * 4] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 0;
        likarray[(i * nSNPs + j) * 4 + 2] = 0;
        likarray[(i * nSNPs + j) * 4 + 3] = 0;
      } else if(SNPs[i*nSNPs + j] == 2) {
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4] = 0;
        likarray[(i * nSNPs + j) * 4 + 2] = 0;
        likarray[(i * nSNPs + j) * 4 + 3] = 0;
      } else if(SNPs[i*nSNPs + j] == 3) {
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 0;
        likarray[(i * nSNPs + j) * 4] = 0;
        likarray[(i * nSNPs + j) * 4 + 3] = 0;
      } else if(SNPs[i*nSNPs + j] == 4) {
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 0;
        likarray[(i * nSNPs + j) * 4] = 0;
        likarray[(i * nSNPs + j) * 4 + 2] = 0;
      } else {
        likarray[(i * nSNPs + j) * 4 + 3] = 1;
        likarray[(i * nSNPs + j) * 4 + 1] = 1;
        likarray[(i * nSNPs + j) * 4 + 2] = 1;
        likarray[(i * nSNPs + j) * 4] = 1;
      }
    }
  }
  for(int i = obs; i < nnodes; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      for(int k = 0; k < 4; ++k) {
        likarray[(i * nSNPs + j) * 4 + k] = 1;
      }
    }
  }
  
  for(int i = 0; i < nlens.size(); ++i) {
    if(nodeparents[i] > 0) {
      nlens[i] = nodetimes[i] - nodetimes[nodeparents[i]-1];
    } else {
      rootnode = i;
    }
  }
  nextnode = 0;
  while(nodeparents[nextnode] - 1 != rootnode) {
    ++nextnode;
  }
  rootnode = nextnode;
  for(int i = 0; i < obs; ++i) {
    routefree[i] = true;
  }
  for(int i = obs; i < 2*obs-1; ++i) {
    routefree[i] = false;
  }
  for(int i = 2*obs - 1; i < routefree.size(); ++i) {
    routefree[i] = true;
  }
  
  for(int i = 0; i < obs; ++i) {
    curnode = i;
    nextnode = nodeparents[curnode] - 1;
    edgelen = nlens[curnode];
    while(routefree[curnode] && nextnode != -1) {
      if(nextnode < 2*obs - 1) {
        for(int j = 0; j < nSNPs; ++j) {
          totprob = likarray[(curnode * nSNPs + j) * 4];
          for(int k = 1; k < 4; ++k) {
            totprob += likarray[(curnode * nSNPs + j) * 4 + k];
          }
          for(int k = 0; k < 4; ++k) {
            likarray[(nextnode * nSNPs + j) * 4 + k] *=
              0.25*totprob + (likarray[(curnode * nSNPs + j) * 4 + k] -
              0.25*totprob)*exp(-mutrate*edgelen);
          }
        }
        curnode = nextnode;
        edgelen = nlens[curnode];
        nextnode = nodeparents[curnode] - 1;
      } else {
        edgelen += nlens[nextnode];
        nextnode = nodeparents[nextnode] - 1;
      }
    }
    routefree[curnode] = true;
  }
  
  for(int j = 0; j < nSNPs; ++j) {
    SNPsums[j] = 0.25*likarray[(rootnode * nSNPs + j) * 4];
    for(int k = 1; k < 4; ++k) {
      SNPsums[j] += 0.25*likarray[(rootnode * nSNPs + j) * 4 + k];
    }
    SNPsums[j] = log(SNPsums[j]) * SNPfreqs[j];
  }
  
  result = SNPsums[0];
  for(int j = 1; j < nSNPs; ++j) {
    result += SNPsums[j];
  }
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
#likseq(t(curstate$d$SNP), curstate$d$SNPfr,
#         curstate$v$nodeparents, curstate$v$nodetimes,.01,200)
*/
