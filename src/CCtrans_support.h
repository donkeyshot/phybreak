#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

std::vector<int> ptroottrans(std::vector<int> pars, int ID);

void makecladearray(std::vector<bool> &ca, const std::vector<int> &pars, 
                    int n, int ss);

void tallyclades(std::vector<int> &pos, std::vector<int> &count,
                 std::vector<int> &refs, const std::vector<bool> &ca,
                 int n, int ss);

void cladetimestats(std::vector<double> &sums, std::vector<double> &ssq,
                    const std::vector<double> &tims, 
                    const std::vector<int> &clIDs, 
                    const std::vector<int> &whclID,
                    const std::vector<bool> &ca, int n, int ss);
