#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
#include "CCtrans_support.h"
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

// std::vector<int> ptroottrans(std::vector<int> , int );
// 
// void makecladearray(std::vector<bool> , const std::vector<int> , 
//                     int , int );
// 
// void tallyclades(std::vector<int> , std::vector<int> ,
//                  std::vector<int> , const std::vector<bool> ,
//                  int , int );
// 
// void cladetimestats(std::vector<double> , std::vector<double> ,
//                     const std::vector<double> , 
//                     const std::vector<int> , 
//                     const std::vector<int> ,
//                     const std::vector<bool> , int , int );


// [[Rcpp::export(name=".CCtranstree")]]
std::vector<double> CCtranstree(const std::vector<int> &pars, 
                             const std::vector<double> &tims,
                             std::vector<int> dims) {

  //For each posterior tree, an obs*obs array of clades is created.
  //Each tree contains 'obs' clades, which are described by vectors
  //of length 'obs'. Each position i in such vector indicates host
  //with ID=i is part of that clade.
  std::vector<bool> cladearray(dims[0] * dims[0] * dims[1]); //dims[0] = obs, dims[1] = no.of.trees
  makecladearray(cladearray, pars, dims[0], dims[1]);
  
  //Clades are identified by their position in cladearray, the first time
  //they appear. The vector 'uniquecladepos' contains these positions,
  //and 'uniquecladecount' records how frequent each clade is observed.
  //Finally, 'whichclade' indicates the cladeID for each clade in cladearray.
  std::vector<int> uniquecladepos;
  std::vector<int> uniquecladecount;
  std::vector<int> whichclade(dims[0] * dims[1]);
  tallyclades(uniquecladepos, uniquecladecount, whichclade,
              cladearray, dims[0], dims[1]);

  //For each unique clade, the sums and sums-of-squares of the infection
  //times of the root hosts are calculated
  std::vector<double> cladetimesums(uniquecladepos.size());
  std::vector<double> cladetimesumsqs(uniquecladepos.size());
  cladetimestats(cladetimesums, cladetimesumsqs, tims,
                 uniquecladepos, whichclade, cladearray, dims[0], dims[1]);

  //For each clade, its log(frequency)
  std::vector<double> cladescores(dims[0] * dims[1]);
  for(int i = 0; i < dims[0]*dims[1]; ++i) {
    cladescores[i] = std::log(uniquecladecount[whichclade[i]]);
  }

  //Determine the posterior tree with highest sum(log(frequency))
  int postsample = 0;
  std::vector<double> treescores(dims[1]);
  for(int i = 0; i < dims[1]; ++i) {
    for(int j = 0; j < dims[0]; ++j) {
      treescores[i] += cladescores[i * dims[0] + j];
    }
    if(treescores[i] > treescores[postsample]) {
      postsample = i;
    }
  }


  //Make vector with results: parents, clade support,
  //tree infection time,
  //mean infection time, SD(infection time)
  std::vector<double> result(5 * dims[0]);
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = pars[dims[0] * postsample + i];
    result[i + dims[0]] = uniquecladecount[whichclade[dims[0] * postsample + i]];
    
    result[i + 2*dims[0]] = 
      cladetimesums[whichclade[dims[0] * postsample + i]] / result[i + dims[0]];
    
    result[i + 3*dims[0]] = std::sqrt((
      cladetimesumsqs[whichclade[dims[0] * postsample + i]] -
        result[i + 2*dims[0]] * result[i + 2*dims[0]] * result[i + dims[0]]
    ) / (result[i + dims[0]] - 1) );
    
    result[i + 4*dims[0]] = 
      tims[dims[0] * postsample + i];
  }
  



  return result;
}


// You can include R code blocks in C++ files processed with sourceCpp
// (useful for testing and development). The R code will be automatically
// run after the compilation.
//

