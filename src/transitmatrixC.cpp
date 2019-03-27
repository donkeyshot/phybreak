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

// [[Rcpp::export]]
NumericVector transitmatrixC(NumericVector x) {
  
  
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
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
timesTwo(42)
*/
