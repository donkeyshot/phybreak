#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
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
                    std::cout << transitmatrix[8*as + 4*cs + 2*gs + ts - 1 +
                      15*(8*af + 4*cf + 2*gf + tf - 1)];
                    std::cout << "\n";
                    
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

// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//
/*** R
transitmatrixC(.001, .1)
***/

