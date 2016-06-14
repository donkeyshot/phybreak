#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
#include "CCtrans_support.h"
using namespace Rcpp;







// [[Rcpp::export(name=".mtcctree")]]
std::vector<double> CCtranstree(const std::vector<int> &pars, 
                                const std::vector<double> &tims,
                                std::vector<int> dims) {
  
  //The ss*n vector with n parents in each of ss trees is transformed into
  //a n_clade*n uniquecladearray, an ss*n vector indicating which n clades
  //are in each of the ss trees, and an n_clade long vector with clade frequencies
  std::vector<bool> uniq_clades;
  std::vector<int> which_clades(dims[0] * dims[1]);
  std::vector<int> clade_freqs;
  int n_clade = 0;
  buildclades(uniq_clades, clade_freqs, which_clades, n_clade,
              pars, dims[0], dims[1]);
  
  //For each clade, its log(frequency)
  std::vector<double> cladescores(dims[0] * dims[1]);
  for(int i = 0; i < dims[0]*dims[1]; ++i) {
    cladescores[i] = std::log(clade_freqs[which_clades[i]]);
  }
  
  //Determine the posterior tree with highest sum(log(frequency))
  int posttree = 0;
  std::vector<double> treescores(dims[1]);
  for(int i = 0; i < dims[1]; ++i) {
    for(int j = 0; j < dims[0]; ++j) {
      treescores[i] += cladescores[i * dims[0] + j];
    }
    if(treescores[i] > treescores[posttree]) {
      posttree = i;
    }
  }
  
  //Determine the posterior clades
  std::vector<int> postclades(dims[0]);
  for(int i = 0; i < dims[0]; ++i) {
    postclades[i] = which_clades[posttree * dims[0] + i];
  }

  //For each posterior clade, the sums and sums-of-squares of the infection
  //times of the root hosts are calculated
  std::vector<double> cladetimesums(dims[0]);
  std::vector<double> cladetimesumsqs(dims[0]);
  cladetimestats(cladetimesums, cladetimesumsqs, tims,
                 postclades, which_clades, uniq_clades, dims[0], dims[1]);
  

  //Make vector with results: parents, clade support,
  //mean infection time, SD(infection time)
  //tree infection time,
  std::vector<double> result(1 + 5 * dims[0]);
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = pars[dims[0] * posttree + i];
    result[i + dims[0]] = clade_freqs[postclades[i]];
    
    result[i + 2*dims[0]] = 
      cladetimesums[i] / result[i + dims[0]];
    
    result[i + 3*dims[0]] = std::sqrt((
      cladetimesumsqs[i] -
        result[i + 2*dims[0]] * result[i + 2*dims[0]] * result[i + dims[0]]
    ) / (result[i + dims[0]] - 1) );
    
    result[i + 4*dims[0]] = 
      tims[dims[0] * posttree + i];
  }
  result[5 * dims[0]] = posttree;
  
  
  
  
  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

