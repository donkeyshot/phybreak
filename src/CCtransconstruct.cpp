#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
#include <utility>
#include <functional>
#include "CCtrans_support.h"
using namespace Rcpp;

// [[Rcpp::plugins("cpp11")]]



  

  




// [[Rcpp::export(name=".CCtranstreeconstruct")]]
std::vector<double> CCtranstreeconstruct(const std::vector<int> &pars, 
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
  

  //The vector 'orderpos' contains the paired 'clade_freqs' and corresponding position in uniq_clades,
  //ordered from most frequent to least frequent
  std::vector<std::pair<int, int> > orderpos(n_clade);
  for(unsigned int i = 0; i < n_clade; ++i) {
    orderpos[i].first = -clade_freqs[i];
    orderpos[i].second = i;
  }
  std::sort(orderpos.begin(), orderpos.end());

  //The vector 'postclades' contains the tree with highest maximum clade support,
  //though not guaranteed so. It is constructed by sequentially testing the next highest
  //supported clade for compatibility with already accepted clades, and accepting
  //if compatible, until the tree is complete. The clades are sorted by root host.
  std::vector<int> postclades(dims[0]);

  selectclades(postclades, uniq_clades, 
                orderpos, dims[0]);
  
  

  
  
  //The 'parents' vector contains the infector for each host, and 'scores' the corresponding score
  std::vector<int> parents(dims[0]);

  findparents(parents, postclades, uniq_clades, dims[0]);
  
  
  //For each posterior clade, the sums and sums-of-squares of the infection
  //times of the root hosts are calculated
  std::vector<double> cladetimesums(dims[0]);
  std::vector<double> cladetimesumsqs(dims[0]);
  cladetimestats(cladetimesums, cladetimesumsqs, tims,
                 postclades, which_clades, uniq_clades, dims[0], dims[1]);
  

  //Make vector with results: parents, clade support,
  //tree infection time,
  //mean infection time, SD(infection time)
  std::vector<double> result(4 * dims[0]);
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = parents[i];
    result[i + dims[0]] = clade_freqs[postclades[i]];
    
    result[i + 2*dims[0]] = 
      cladetimesums[i] / result[i + dims[0]];
    
    result[i + 3*dims[0]] = std::sqrt((
      cladetimesumsqs[i] -
        result[i + 2*dims[0]] * result[i + 2*dims[0]] * result[i + dims[0]]
    ) / (result[i + dims[0]] - 1) );
    
  }
  
//   std::vector<double> result2(20);
//   for(int i = 0; i < 20; ++i) {
//     result2[i] = parents[i];
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

