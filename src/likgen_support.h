#include <Rcpp.h>
#include <vector>
#include <string>
#include <iterator>
#include <cmath>
using namespace Rcpp;


NumericVector obsarrayC(IntegerVector SNPs, int Nsamples);

NumericVector obshostarrayC(NumericVector obssamplearray, 
                            IntegerVector samplehosts, int Nhosts); 


double facultyC(double x);

double chooseC(double a, double b);



NumericVector transitmatrixC(double mutrate, double bnsize);