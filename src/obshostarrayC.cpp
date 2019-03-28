// #include <Rcpp.h>
// using namespace Rcpp;
// 
// // This is a simple example of exporting a C++ function to R. You can
// // source this function into an R session using the Rcpp::sourceCpp 
// // function (or via the Source button on the editor toolbar). Learn
// // more about Rcpp at:
// //
// //   http://www.rcpp.org/
// //   http://adv-r.had.co.nz/Rcpp.html
// //   http://gallery.rcpp.org/
// //
// 
// // [[Rcpp::export]]
// NumericVector obshostarrayC(NumericVector obssamplearray, 
//                            IntegerVector samplehosts, int Nhosts) {
//   int Nsamples = samplehosts.size();
//   int nSNPs = obssamplearray.size() / Nsamples / 15;
//   NumericVector obshostarray(Nhosts * nSNPs * 15);
//   
//   for(int i = 0; i < obshostarray.size(); ++i) {
//     obshostarray[i] = 1;
//   }
//   for(int sample = 0; sample < Nsamples; ++sample) {
//     for(int nuc = 0; nuc < nSNPs * 15; ++nuc) {
//       obshostarray[(samplehosts[sample] - 1) * nSNPs * 15 + nuc] *=
//         obssamplearray[sample * nSNPs * 15 + nuc];
//     }
//   }
//   
//   return obshostarray;
// }
// 
// 
// // You can include R code blocks in C++ files processed with sourceCpp
// // (useful for testing and development). The R code will be automatically 
// // run after the compilation.
// //
// 
// /*** R
// */
