#include <Rcpp.h>
#include <vector>
// #include <iterator>
// #include <cmath>
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

// [[Rcpp::export(name=".ptr")]]
std::vector<int> ptr(IntegerVector pars, int ID) {
  std::vector<int> ans;
  ans.push_back(ID);
  int nextID = pars[ID - 1];
  while(nextID > 0) {
    ans.push_back(nextID);
    nextID = pars[nextID - 1];
  }
  return(ans);
}




