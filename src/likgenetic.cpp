#include <Rcpp.h>
#include <vector>
#include <string>
using namespace Rcpp;

// Calculation of log-likelihood of sequence data
// given a phylo-transmission tree

// [[Rcpp::export(name=".likgenetic")]]
double likgenetic(IntegerVector SNPs, IntegerVector SNPfreqs,
              IntegerVector infectors, NumericVector samplehosts,
              double mutrate, double bnsize) {
  int Nsamples = samplehosts.size();
  int Nhosts = infectors.size();
  int nnodes = Nsamples + Nhosts;
  int nSNPs = SNPs.size()/Nsamples;
  NumericVector obsarray(Nsamples * nSNPs * 16);
  NumericVector likarray(Nhosts * nSNPs * 16);
  
  NumericVector changearray(4 * 4 * 4);
  NumericVector SNPsums(nSNPs);
  NumericVector routeblock(Nhosts);
  int curnode;
  int nextnode;
  int indexhost = 0;
  double edgelen;
  double totprob;
  double result;
  
  // observation probabilities
  obsarray = obsarrayC(SNPs, Nsamples);
  
  for(int i = 0; i < Nhosts; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      for(int k = 0; k < 16; ++k) {
        likarray[(i * nSNPs + j) * 16 + k] = 1;
      }
    }
  }
  
  // transition matrix
  for(int i = 0; i < 4; ++i) {
    for(int j = 0; j < i + 1; ++j) {
      for(int k = j; k < j + 4 - i; ++k) {
        changearray[i * 16 + j * 4 + k] = 
          (1 - bnsize) * ((j + 1) / (i + 1)) * 
          mutrate ^ k * (1 - mutrate) ^ (3 - k);
        for(int bn = 1; bn < j + 1) {
          changearray[i * 16 + j * 4 + k] +=
            bnsize ^ bn * ()
        }
        cstart <- c(rep(1, i), rep(0, statecount - i))
        cfinish <- c(rep(1, j), rep(0, statecount - k), rep(1, k - j))
        changearray[i * 16 + j * 4 + k] <-
          pr_change(cstart, cfinish, mu, bn)
      }
    }
  }
  sum(
    sapply(1:sum(c_start),
           function(x) {
      bn ^ (x - 1) *
        (1 - bn) ^ (sum(c_start) > x) *
        (choose(sum(c_start * c_finish), x) /
          choose(sum(c_start), x)) *
            mu ^ (sum(c_finish) - x) *
            (1 - mu) ^ (length(c_start) - sum(c_finish))
    })
  )
    
  
  
 // nextnode = 0;
 // while(nodeparents[nextnode] - 1 != rootnode) {
 //   ++nextnode;
 // }
 // rootnode = nextnode;
  for(int i = 0; i < Nhosts; ++i) {
    routeblock[i] = 0;
  }
  for(int i = 0; i < Nhosts; ++i) {
    if(infectors[i] == 0) {
      indexhost = i;
    } else {
      routeblock[infectors[i] - 1] += 1;
    }
  }

  
  for(int i = 0; i < Nsamples; ++i) {
      for(int j = 0; j < nSNPs; ++j) {
        totprob = likarray[(curnode * nSNPs + j) * 4];
        for(int k = 1; k < 4; ++k) {
          totprob += likarray[(curnode * nSNPs + j) * 4 + k];
        }
        for(int k = 0; k < 4; ++k) {
          likarray[(nextnode * nSNPs + j) * 4 + k] *=
            0.25 * totprob + (likarray[(curnode * nSNPs + j) * 4 + k] -
            0.25 * totprob) * exp(-mutrate * edgelen);
        }
      }
      curnode = nextnode;
      edgelen = nlens[curnode];
      nextnode = nodeparents[curnode] - 1;
      //      } else {
      //        edgelen += nlens[nextnode];
      //        nextnode = nodeparents[nextnode] - 1;
      //      }
    }
    routefree[curnode] = true;
  }
  
  
  for(int i = 0; i < Nsamples; ++i) {
    curnode = i;
    nextnode = nodeparents[curnode] - 1;
    edgelen = nlens[curnode];
    while(routefree[curnode] && nextnode != -1) {
//      if(nextnode < 2*Nsamples - 1) {
        for(int j = 0; j < nSNPs; ++j) {
          totprob = likarray[(curnode * nSNPs + j) * 4];
          for(int k = 1; k < 4; ++k) {
            totprob += likarray[(curnode * nSNPs + j) * 4 + k];
          }
          for(int k = 0; k < 4; ++k) {
            likarray[(nextnode * nSNPs + j) * 4 + k] *=
              0.25 * totprob + (likarray[(curnode * nSNPs + j) * 4 + k] -
              0.25 * totprob) * exp(-mutrate * edgelen);
          }
        }
        curnode = nextnode;
        edgelen = nlens[curnode];
        nextnode = nodeparents[curnode] - 1;
//      } else {
//        edgelen += nlens[nextnode];
//        nextnode = nodeparents[nextnode] - 1;
//      }
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
