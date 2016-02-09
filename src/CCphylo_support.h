#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

std::vector<int> ptrootphylo(std::vector<int> pars, int ID);

void makephylocladearray(std::vector<bool> &ca, const std::vector<int> &pars, 
                         int n, int ss);

void tallyphyloclades(std::vector<int> &pos, std::vector<int> &count,
                      std::vector<int> &refs, const std::vector<bool> &ca,
                      int n, int ss);

void buildphyloclades(std::vector<bool> &uc, std::vector<int> &cf, std::vector<int> &wc, 
                      int &nc, const std::vector<int> &pars, int n, int ss);

void phylocladetimestats(std::vector<double> &sums, std::vector<double> &ssq,
                         const std::vector<double> &tims, 
                         const std::vector<int> &pc, 
                         const std::vector<int> &wc,
                         const std::vector<bool> &uc, int n, int ss);


//pc = postclade, will contain n positions of clades in 'uc' that are accepted 
//uc = array with all unique clades, with each n bool to indicate the hosts
//op = ordered pairs of clade indicators (positions in 'uc') and their frequencies (negative) among all posterior trees
//n = outbreak size
void selectphyloclades(std::vector<int> &pc, const std::vector<bool> &uc, 
                       const std::vector<std::pair<int, int> > &op, int n);


//Function to find the node parents given the bestcladetree parents, scores, bestcladetree, besttreescores, cladeparents
void findnodeparents(std::vector<int> &pars, const std::vector<int> &pc, 
                     const std::vector<bool> &uc, int n);
