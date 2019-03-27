#include <Rcpp.h>
#include <vector>
#include <string>
#include <iterator>
#include <cmath>
using namespace Rcpp;


// [[Rcpp::export]]
NumericVector obsarrayC(IntegerVector SNPs, int Nsamples);

// [[Rcpp::export]]
NumericVector obshostarrayC(NumericVector obssamplearray, 
                            IntegerVector samplehosts, int Nhosts); 


// [[Rcpp::export]]
double facultyC(double x);

// [[Rcpp::export]]
double chooseC(double a, double b);



// [[Rcpp::export]]
NumericVector transitmatrixC(double mutrate, double bnsize);