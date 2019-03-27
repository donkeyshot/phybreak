#include <Rcpp.h>
using namespace Rcpp;

// This is a simple example of exporting a C++ function to R. You can
// source this function into an R session using the Rcpp::sourceCpp 
// function (or via the Source button on the editor toolbar). Learn
// more about Rcpp at:
//
//   http://www.rcpp.org/
//   http://adv-r.had.co.nz/Rcpp.html
//   http://gallery.rcpp.org/
//

// [[Rcpp::export(name=.obsarray)]]
NumericVector obsarrayC(IntegerVector SNPs, int Nsamples) {
  int nSNPs = SNPs.size()/Nsamples;
  int apres;
  int cpres;
  int gpres;
  int tpres;
  NumericVector obsarray(Nsamples * nSNPs * 16);
  
  for(int i = 0; i < Nsamples; ++i) {
    for(int j = 0; j < nSNPs; ++j) {
      obsarray[(i * nSNPs + j) * 16] = 0;
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
      for(int nuca = 0; nuca < 2; ++nuca) {
        for(int nucc = 0; nucc < 2; ++nucc) {
          for(int nucg = 0; nucg < 2; ++nucg) {
            for(int nuct = 0; nuct < 2; ++nuct) {
              if(nuca + nucc + nucg + nuct > 0) {
                obsarray[(i * nSNPs + j) * 16 +
                  nuca * 8 + nucc * 4 + nucg * 2 + nuct] = 
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


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

