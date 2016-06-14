// #include <Rcpp.h>
// #include <vector>
// #include <iterator>
// #include <cmath>
// #include "CCphylo_support.h"
// using namespace Rcpp;
// 
// 
// 
// // [[Rcpp::export(name=".CCphylotree")]]
// std::vector<double> CCphylotree(const std::vector<int> &pars, 
//                                 const std::vector<double> &tims,
//                                std::vector<int> dims) {
// 
//   
//   //The ss*n vector with n parents in each of ss trees is transformed into
//   //a n_clade*n uniquecladearray, an ss*(n-1) vector indicating which n-1 clades
//   //are in each of the ss trees, and an n_clade long vector with clade frequencies
//   std::vector<bool> uniq_clades;
//   std::vector<int> which_clades((dims[0] - 1) * dims[1]);
//   std::vector<int> clade_freqs;
//   int n_clade = 0;
//   buildphyloclades(uniq_clades, clade_freqs, which_clades, n_clade,
//               pars, dims[0], dims[1]);
//   
//   
//   //For each clade, its log(frequency)
//   std::vector<double> cladescores((dims[0] - 1) * dims[1]);
//   for(int i = 0; i < (dims[0]-1)*dims[1]; ++i) {
//     cladescores[i] = std::log(clade_freqs[which_clades[i]]);
//   }
// 
//   //Determine the posterior tree with highest sum(log(frequency))
//   int posttree = 0;
//   std::vector<double> treescores(dims[1]);
//   for(int i = 0; i < dims[1]; ++i) {
//     for(int j = 0; j < dims[0] - 1; ++j) {
//       treescores[i] += cladescores[i * (dims[0] - 1) + j];
//     }
//     if(treescores[i] > treescores[posttree]) {
//       posttree = i;
//     }
//   }
//   
//   //Determine the posterior clades
//   std::vector<int> postclades(dims[0] - 1);
//   for(int i = 0; i < dims[0] - 1; ++i) {
//     postclades[i] = which_clades[posttree * (dims[0] - 1) + i];
//   }
//   
//   //For each posterior clade, the sums and sums-of-squares of the infection
//   //times of the root hosts are calculated
//   std::vector<double> cladetimesums(dims[0] - 1);
//   std::vector<double> cladetimesumsqs(dims[0] - 1);
//   phylocladetimestats(cladetimesums, cladetimesumsqs, tims,
//                  postclades, which_clades, uniq_clades, dims[0], dims[1]);
//   
// 
//   //Make vector with results: parents, clade support,
//   //mean node time, SD(node time)
//   //tree infection time,
//   std::vector<double> result(5 * (2*dims[0] - 1));
//   for(int i = 0; i < dims[0]; ++i) {
//     result[i] = pars[(2*dims[0]-1) * posttree + i];
//     result[i + 2*dims[0] - 1] = dims[1];
//     result[i + 2* (2*dims[0] - 1)] = 0;
//     result[i + 3* (2*dims[0] - 1)] = 0;
//     result[i + 4* (2*dims[0] - 1)] = 0;
//   }
//   for(int i = 0; i < dims[0] - 1; ++i) {
//     result[i + dims[0]] = pars[(2*dims[0]-1) * posttree + dims[0] + i];
//     result[i + dims[0] + 2*dims[0] - 1] = clade_freqs[postclades[i]];
//     result[i + dims[0] + 2*(2*dims[0] - 1)] = (
//       cladetimesums[i] / result[i + dims[0] + 2*dims[0] - 1]);
//     result[i + dims[0] + 3*(2*dims[0] - 1)] = std::sqrt((
//       cladetimesumsqs[i] -
//         result[i + dims[0] + 2*(2*dims[0] - 1)] * 
//         result[i + dims[0] + 2*(2*dims[0] - 1)] * 
//         result[i + dims[0] + 2*dims[0] - 1]
//     ) / (result[i + dims[0] + 2*dims[0] - 1] - 1) );
//     result[i + dims[0] + 4*(2*dims[0] - 1)] = tims[(dims[0]-1) * posttree + i];
//     
//   }
//       
//       
// 
// 
//   return result;
// }
// 
// 
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically
// // run after the compilation.
// //
// 
