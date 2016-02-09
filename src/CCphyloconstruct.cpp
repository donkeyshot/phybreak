#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
#include "CCphylo_support.h"
using namespace Rcpp;




// [[Rcpp::export(name=".CCphylotreeconstruct")]]
std::vector<double> CCphyloconstruct(const std::vector<int> &pars, 
                                  const std::vector<double> &tims,
                                  std::vector<int> dims) {
  
  
  //The ss*n vector with n parents in each of ss trees is transformed into
  //a n_clade*n uniquecladearray, an ss*(n-1) vector indicating which n-1 clades
  //are in each of the ss trees, and an n_clade long vector with clade frequencies
  std::vector<bool> uniq_clades;
  std::vector<int> which_clades((dims[0] - 1) * dims[1]);
  std::vector<int> clade_freqs;
  int n_clade = 0;
  buildphyloclades(uniq_clades, clade_freqs, which_clades, n_clade,
                   pars, dims[0], dims[1]);
  
  //The vector 'orderpos' contains the paired 'clade_freqs' and corresponding position in uniq_clades,
  //ordered from most frequent to least frequent
  std::vector<std::pair<int, int> > orderpos(n_clade);
  for(int i = 0; i < n_clade; ++i) {
    orderpos[i].first = -clade_freqs[i];
    orderpos[i].second = i;
  }
  std::sort(orderpos.begin(), orderpos.end());
  
  
  //The vector 'postclades' contains the tree with highest maximum clade support,
  //though not guaranteed so. It is constructed by sequentially testing the next highest
  //supported clade for compatibility with already accepted clades, and accepting
  //if compatible, until the tree is complete. 
  std::vector<int> postclades(dims[0]);
  
  selectphyloclades(postclades, uniq_clades, 
               orderpos, dims[0]);
  
  //The 'parents' vector contains the parent for each node, and 'scores' the corresponding score
  std::vector<int> parents(2*dims[0] - 1);
  
  findnodeparents(parents, postclades, uniq_clades, dims[0]);
  
  //For each posterior clade, the sums and sums-of-squares of the infection
  //times of the root hosts are calculated
  std::vector<double> cladetimesums(dims[0] - 1);
  std::vector<double> cladetimesumsqs(dims[0] - 1);
  phylocladetimestats(cladetimesums, cladetimesumsqs, tims,
                      postclades, which_clades, uniq_clades, dims[0], dims[1]);
  

  //Make vector with results: parents, clade support,
  //mean node time, SD(node time)
  //tree infection time,
  std::vector<double> result(4 * (2*dims[0] - 1));
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = parents[i];
    result[i + 2*dims[0] - 1] = dims[1];
    result[i + 2* (2*dims[0] - 1)] = 0;
    result[i + 3* (2*dims[0] - 1)] = 0;
  }
  for(int i = 0; i < dims[0] - 1; ++i) {
    result[i + dims[0]] = parents[i + dims[0]];
    result[i + dims[0] + 2*dims[0] - 1] = clade_freqs[postclades[i]];
    result[i + dims[0] + 2*(2*dims[0] - 1)] = (
      cladetimesums[i] / result[i + dims[0] + 2*dims[0] - 1]);
    result[i + dims[0] + 3*(2*dims[0] - 1)] = std::sqrt((
      cladetimesumsqs[i] -
        result[i + dims[0] + 2*(2*dims[0] - 1)] * 
        result[i + dims[0] + 2*(2*dims[0] - 1)] * 
        result[i + dims[0] + 2*dims[0] - 1]
    ) / (result[i + dims[0] + 2*dims[0] - 1] - 1) );

  }
  //
//   std::vector<int> result2(2450);
//   for(int i = 0; i < 2450; ++i) {
//     result2[i] = bestcladetree[i];
//   }
  //     std::vector<int> result2(200);
  //     for(int i = 0; i < 100; ++i) {
  //       result2[2*i] = orderpos[i].first;
  //       result2[2*i+1] = orderpos[i].second;
  //     }



  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

