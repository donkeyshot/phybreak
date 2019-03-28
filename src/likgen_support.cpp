#include <Rcpp.h>
#include <vector>
#include <string>
#include <iterator>
#include <cmath>
#include "likgen_support.h"
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector obsarrayC(IntegerVector SNPs, int Nsamples) {
  int nSNPs = SNPs.size()/Nsamples;
  double apres;
  double cpres;
  double gpres;
  double tpres;
  NumericVector obsarray(Nsamples * nSNPs * 15);

  for(int i = 0; i < Nsamples; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      obsarray[(i * nSNPs + j) * 15] = 0;
      apres = 0;
      cpres = 0;
      gpres = 0;
      tpres = 0;
      if(SNPs[i * nSNPs + j] == 1) {
        apres = 1;
      }
      if(SNPs[i * nSNPs + j] == 2) {
        cpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 3) {
        gpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 4) {
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 5) {
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 6) {
        apres = 1;
        cpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 7) {
        apres = 1;
        gpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 8) {
        apres = 1;
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 9) {
        cpres = 1;
        gpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 10) {
        cpres = 1;
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 11) {
        gpres = 1;
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 12) {
        apres = 1;
        cpres = 1;
        gpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 13) {
        apres = 1;
        cpres = 1;
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 14) {
        apres = 1;
        gpres = 1;
        tpres = 1;
      }
      if(SNPs[i * nSNPs + j] == 15) {
        cpres = 1;
        gpres = 1;
        tpres = 1;
      }
      for(double nuca = 0; nuca < 2; ++nuca) {
        for(double nucc = 0; nucc < 2; ++nucc) {
          for(double nucg = 0; nucg < 2; ++nucg) {
            for(double nuct = 0; nuct < 2; ++nuct) {
              if(nuca + nucc + nucg + nuct > 0) {
                obsarray[(i * nSNPs + j) * 15 +
                  nuca * 8 + nucc * 4 + nucg * 2 + nuct - 1] =
                  (nuca * apres + nucc * cpres + nucg * gpres +
                  nuct * tpres) / (nuca + nucc + nucg + nuct);
              }
            }
          }
        }
      }
    }
  }

  return obsarray;
}

// [[Rcpp::export]]
NumericVector obshostarrayC(NumericVector obssamplearray,
                            IntegerVector samplehosts, int Nhosts) {
  int Nsamples = samplehosts.size();
  int nSNPs = obssamplearray.size() / Nsamples / 15;
  NumericVector obshostarray(Nhosts * nSNPs * 15);

  for(int i = 0; i < obshostarray.size(); ++i) {
    obshostarray[i] = 1;
  }
  for(int sample = 0; sample < Nsamples; ++sample) {
    for(int nuc = 0; nuc < nSNPs * 15; ++nuc) {
      obshostarray[(samplehosts[sample] - 1) * nSNPs * 15 + nuc] *=
        obssamplearray[sample * nSNPs * 15 + nuc];
    }
  }

  return obshostarray;
}


// [[Rcpp::export]]
double facultyC(double x) {
  double res;

  if(x == 0) {
    res = 1;
  } else {
    res = x * facultyC(x - 1);
  }

  return res;
}

// [[Rcpp::export]]
double chooseC(double a, double b) {
  double res;

  if(b > a) {
    res = 0;
  } else {
    res = facultyC(a) / (facultyC(b) * (facultyC(a - b)));
  }

  return res;
}



// [[Rcpp::export]]
NumericVector transitmatrixC(double mutrate, double bnsize) {
  NumericVector transitmatrix(15 * 15);
  double bottleneck;
  double nstart;
  double nfinish;
  double nsame;

  for(int as = 0; as < 2; ++as) {
    for(int cs = 0; cs < 2; ++cs) {
      for(int gs = 0; gs < 2; ++gs) {
        for(int ts = 0; ts < 2; ++ts) {
          for(int af = 0; af < 2; ++af) {
            for(int cf = 0; cf < 2; ++cf) {
              for(int gf = 0; gf < 2; ++gf) {
                for(int tf = 0; tf < 2; ++tf) {
                  nstart = as + cs + gs + ts;
                  nfinish = af + cf + gf + tf;
                  nsame = as*af + cs*cf + gs*gf + ts*tf;
                  if(nsame > 0) {
                    if(nstart > 1) {
                      bottleneck = 1 - bnsize;
                    } else {
                      bottleneck = 1;
                    }
                    transitmatrix[8*as + 4*cs + 2*gs + ts - 1 +
                      15*(8*af + 4*cf + 2*gf + tf - 1)] +=
                      bottleneck * (nsame/nstart) *
                      std::pow(mutrate, (nfinish - 1)) * std::pow((1-mutrate),(4-nfinish));
                    for(int i = 1; i < nsame; ++i) {
                      if(nstart > i + 1) {
                        bottleneck = 1 - bnsize;
                      } else {
                        bottleneck = 1;
                      }
                      transitmatrix[8*as + 4*cs + 2*gs + ts - 1 +
                        15*(8*af + 4*cf + 2*gf + tf - 1)] +=
                        std::pow(bnsize,i) * bottleneck *
                        (chooseC(nsame,i+1)/chooseC(nstart,i+1)) *
                        std::pow(mutrate,(nfinish-i-1)) * std::pow((1-mutrate),(4-nfinish));

                    }
                  }
                  // else {
                  //   if(nstart * nfinish > 0) {
                  //     transitmatrix[8*as + 4*cs + 2*gs + ts - 1 +
                  //       15*(8*af + 4*cf + 2*gf + tf - 1)] = 0;
                  //   }
                  // }
                }
              }
            }
          }
        }
      }
    }
  }
  return transitmatrix;
}
