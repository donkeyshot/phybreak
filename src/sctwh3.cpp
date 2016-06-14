#include <Rcpp.h>
#include <vector>
#include <functional>
#include <algorithm>
using namespace Rcpp;

// Sample coalescence times


// [[Rcpp::export(name=".sctwh3")]]
NumericVector sctwh3(NumericVector tle) { //in reverse order
  int nleafs = tle.size();
  std::vector<long double> tnodetrans;
  std::vector<long double> tintervals;
  std::vector<long double> cumrates;
  long double curleaftime;
  unsigned int nthnode;
  long double newcumrate;
  long double randnr;
  long double newnodetime;
  
  randnr = -log(runif(1,0,1)[0]);
  
  tnodetrans.push_back(tle[1]-randnr);
  
  for(int i = 2; i < nleafs; ++i) {
    curleaftime = tle[i];
    if(curleaftime < tnodetrans[0]) {
      cumrates.push_back(curleaftime);
      nthnode = 0;
    } else {
      cumrates.push_back(tnodetrans[0]);
      
      nthnode = 1;
      while(nthnode < tnodetrans.size() && curleaftime > tnodetrans[nthnode]) {
        newcumrate = (tnodetrans[nthnode] - tnodetrans[nthnode - 1]) * (nthnode + 1);
        cumrates.push_back(newcumrate);
        ++nthnode;
      }
      newcumrate = (curleaftime - tnodetrans[nthnode - 1]) * (nthnode + 1);
      cumrates.push_back(newcumrate);
    }
    
    randnr = -log(runif(1,0,1)[0]);
    tintervals = tnodetrans;
    tintervals.resize(nthnode);
    tintervals.push_back(curleaftime);
    
    while(nthnode > 0 && randnr > cumrates[nthnode]) {
      randnr -= cumrates[nthnode];
      --nthnode;
    }
    
    newnodetime = tintervals[nthnode] - randnr / (nthnode + 1);
    tnodetrans.push_back(newnodetime);
    sort(tnodetrans.begin(),tnodetrans.end());
    
    cumrates.resize(0);
    
  }
  
  NumericVector res(tnodetrans.size());
  for(unsigned int i = 0; i < tnodetrans.size(); ++i) {
    res[i] = tnodetrans[i];
  }
  return res;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

/*** R
.sctwh3(c(4,3,2,1))
*/
