#include <Rcpp.h>
#include <vector>
#include <string>
#include "likgen_support.h"
using namespace Rcpp;

// Calculation of log-likelihood of sequence data
// given a phylo-transmission tree

// [[Rcpp::export]]
double likgeneticC(IntegerVector SNPs, IntegerVector SNPfreqs,
              IntegerVector infectors, IntegerVector samplehosts,
              double mutrate, double bnsize) {
  int Nsamples = samplehosts.size();
  int Nhosts = infectors.size();
  int nSNPs = SNPs.size()/Nsamples;
  NumericVector obssamplearray(Nsamples * nSNPs * 15);
  NumericVector probleaves(Nhosts * nSNPs * 15);
  NumericVector transitmatrix(15 * 15);
  NumericVector likarray(Nhosts * nSNPs * 15);
  mutrate = 1 - exp(-mutrate);
  bnsize = 1 - exp(-bnsize);
  
  NumericVector SNPsums(nSNPs);
  NumericVector leavestogo(Nhosts);
  NumericVector priorprobs(15);
  int curhost;
  int nexthost;
  int indexhost = 0;
  double result;
  
  
  // conditional observation probabilities per sample
  obssamplearray = obsarrayC(SNPs, Nsamples);
  
  // conditional observation probabilities per host
  probleaves = obshostarrayC(obssamplearray, samplehosts, Nhosts);
  
  // transition matrix (from row to column)
  transitmatrix = transitmatrixC(mutrate, bnsize);
   

 // nextnode = 0;
 // while(nodeparents[nextnode] - 1 != rootnode) {
 //   ++nextnode;
 // }
 // rootnode = nextnode;
  for(int i = 0; i < Nhosts; ++i) {
    leavestogo[i] = 0;
  }
  for(int i = 0; i < Nhosts; ++i) {
    if(infectors[i] == 0) {
      indexhost = i;
    } else {
      leavestogo[infectors[i] - 1] += 1;
    }
  }

  
  for(int i = 0; i < Nhosts; ++i) {
    curhost = i;
    nexthost = infectors[i] - 1;
    while(curhost != -1 && leavestogo[curhost] == 0) {
      for(int j = 0; j < nSNPs; ++j) {
        for(int k = 0; k < 15; ++k) {
          for(int l = 0; l < 15; ++l) {
            likarray[curhost * nSNPs * 15 + j * 15 + k] +=
              transitmatrix[15 * l + k] * 
              probleaves[curhost * nSNPs * 15 + j * 15 + l];
          }
        }
      }
      if(nexthost != -1) {
        for(int nuc = 0; nuc < nSNPs * 15; ++nuc) {
          probleaves[nexthost * nSNPs * 15 + nuc] *=
            likarray[curhost * nSNPs * 15 + nuc];
        }
      }
      leavestogo[curhost]--;
      curhost = nexthost;
      nexthost = infectors[curhost] - 1;
      leavestogo[curhost]--;
    }
  }

  for(int k = 0; k < 15; ++k) {
    priorprobs[k] =
      (transitmatrix[k * 15] + transitmatrix[k * 15 + 1] +
      transitmatrix[k * 15 + 3] + transitmatrix[k * 15 + 7])/4;
  }

  for(int j = 0; j < nSNPs; ++j) {
    for(int k = 0; k < 15; ++k) {
      SNPsums[j] += priorprobs[k] * 
        likarray[(indexhost * nSNPs * 15 + j * 15 + k)];
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
