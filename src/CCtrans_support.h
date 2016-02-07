#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

std::vector<int> ptroottrans(std::vector<int> pars, int ID) ; 


void makecladearray(std::vector<bool> &ca, const std::vector<int> &pars, 
                    int n, int ss) ;


void tallyclades(std::vector<int> &pos, std::vector<int> &count,
                 std::vector<int> &refs, const std::vector<bool> &ca,
                 int n, int ss) ;

void buildclades(std::vector<bool> &uc, std::vector<int> &cf, std::vector<int> &wc, 
                 int &nc, const std::vector<int> &pars, int n, int ss) ;


void cladetimestats(std::vector<double> &sums, std::vector<double> &ssq,
                    const std::vector<double> &tims, 
                    const std::vector<int> &pc, 
                    const std::vector<int> &wc,
                    const std::vector<bool> &uc, int n, int ss);


void selectclades(std::vector<int> &pc, const std::vector<bool> &uc, 
                  const std::vector<std::pair<int, int> > &op, int n) ;


void findparents(std::vector<int> &pars, const std::vector<int> &pc, 
                 const std::vector<bool> &uc, int n);

