#include <Rcpp.h>
#include <vector>
#include <iterator>
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

// std::vector<int> ptrootphylo(std::vector<int> pars, int ID) {
//   std::vector<int> ans;
//   int nextID = pars[ID - 1];
//   while(nextID > 0) {
//     ans.push_back(nextID);
//     nextID = pars[nextID - 1];
//   }
//   return(ans);
// }

// std::vector<int> ptroottrans(std::vector<int> pars, int ID) {
//   std::vector<int> ans;
//   ans.push_back(ID);
//   int nextID = pars[ID - 1];
//   while(nextID > 0) {
//     ans.push_back(nextID);
//     nextID = pars[nextID - 1];
//   }
//   return(ans);
// }


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically 
// run after the compilation.
//

/*** R
# timesTwo(42)
*/
