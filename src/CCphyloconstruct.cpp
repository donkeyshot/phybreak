#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
using namespace Rcpp;

// Calculation of the path to the root in a tree,
// given a parent vector and node/host ID

std::vector<int> ptroot1(std::vector<int> pars, int ID) {
  std::vector<int> ans;
  int nextID = pars[ID - 1];
  while(nextID > 0) {
    ans.push_back(nextID);
    nextID = pars[nextID - 1];
  }
  return(ans);
}




// [[Rcpp::export(name=".CCphylotreeconstruct")]]
std::vector<int> CCphyloconstruct(std::vector<int> pars, std::vector<int> dims) {
  std::vector<bool> cladearray(dims[0] * (dims[0] - 1) * dims[1]); //dims[0] = obs, dims[1] = samplesize
  std::vector<int> parset(2 * dims[0] - 1);
  std::vector<int> pos;

  for(int i = 0; i < dims[1]; ++i) {
    for(int ii = 0; ii< 2*dims[0] - 1; ++ii) {
      parset[ii] = pars[i*(2*dims[0]-1) + ii];
    }
    for(int j = 0; j < dims[0]; ++j) {
      pos = ptroot1(parset, j + 1);
      for(unsigned int k = 0; k < pos.size(); ++k) {
        cladearray[i*dims[0]*(dims[0]-1) + (pos[k]-1-dims[0])*dims[0] + j] = true;
      }
    }
  }

  std::vector<int> uniquecladepos;
  std::vector<int> uniquecladecount;
  std::vector<int> whichclade((dims[0]-1) * dims[1]);
  bool cladeequal = false;
  bool nodeequal = true;
  unsigned int nextclade = 0;
  int nextnode = 0;
  uniquecladepos.push_back(0);
  uniquecladecount.push_back(1);
  for(int i = 1; i < (dims[0] - 1) * dims[1]; ++i) {
    while(!cladeequal && nextclade < uniquecladepos.size()) {
      while(nodeequal && nextnode < dims[0]) {
        if(cladearray[i*dims[0] + nextnode] !=
           cladearray[uniquecladepos[nextclade]*dims[0] + nextnode]) {
          nodeequal = false;
        }
        ++nextnode;
      }
      cladeequal = nodeequal;
      nodeequal = true;
      ++nextclade;
      nextnode = 0;
    }
    if(!cladeequal) {
      uniquecladepos.push_back(i);
      uniquecladecount.push_back(1);
      whichclade[i] = uniquecladepos.size() - 1;
    } else {
      ++uniquecladecount[nextclade - 1];
      whichclade[i] = nextclade - 1;
    }
    cladeequal = false;
    nodeequal = true;
    nextclade = 0;
    nextnode = 0;
  }

  std::vector<std::pair<int, int> > orderpos(uniquecladepos.size());
  for(unsigned int i = 0; i < orderpos.size(); ++i) {
    orderpos[i].first = -uniquecladecount[i];
    orderpos[i].second = uniquecladepos[i];
  }
  std::sort(orderpos.begin(), orderpos.end());


  std::vector<bool> bestcladetree(dims[0] * (dims[0] - 1));
  std::vector<int> besttreescores(dims[0]);
  for(int i = 0; i < dims[0]; ++i) {
    bestcladetree[i] = cladearray[orderpos[0].second * dims[0] + i];
    bestcladetree[i + dims[0]] =
      cladearray[orderpos[1].second * dims[0] + i];
  }
  besttreescores[0] = -orderpos[0].first;
  besttreescores[1] = -orderpos[1].first;

  nextclade = 2;
  int cladecount = 2;
  bool cladeaccept;
  bool overlap, newmore, bestmore;
  int newclade = 0;


  while(cladecount < dims[0] - 1) {
    newclade = orderpos[nextclade].second * dims[0];
    cladeaccept = true;
    for(int i = 0; i < cladecount * (dims[0] - 1) && cladeaccept; i += dims[0]) {
      overlap = false;
      newmore = false;
      bestmore = false;
      for(int j = 0; j < dims[0]; ++j) {
        if(cladearray[newclade + j] && bestcladetree[i + j]) {overlap = true;}
        if(cladearray[newclade + j] && !bestcladetree[i + j]) {newmore = true;}
        if(!cladearray[newclade + j] && bestcladetree[i + j]) {bestmore = true;}
      }
      if(overlap && newmore && bestmore) {
        cladeaccept = false;
      } else {cladeaccept = true;}
    }
    if(cladeaccept) {
      for(int j = 0; j < dims[0]; ++j) {
        bestcladetree[cladecount * dims[0] + j] = cladearray[newclade + j];
      }
      besttreescores[cladecount] = -orderpos[nextclade].first;
      ++cladecount;
    }
    ++nextclade;
  }

  std::vector<int> cladesizes(dims[0] + 1);
  std::vector<int> parents(2*dims[0]-1);
  int smallest = 0;
  int secsmallest = 0;
  bool cladeinclade = false;
  for(int i = 0; i < dims[0] - 1; ++i) {
    for(int j = 0; j < dims[0]; ++j) {
      cladesizes[i] += bestcladetree[i * dims[0] + j];
    }
  }
  cladesizes[dims[0] - 1] = dims[0] + 1;
  cladesizes[dims[0]] = dims[0] + 2;
  for(int i = 0; i < dims[0]; ++i) {
    smallest = dims[0] - 1;
    for(int j = 0; j < dims[0] - 1; ++j) {
      if(bestcladetree[j * dims[0] + i]) {
        if(cladesizes[j] < cladesizes[smallest]) {
          smallest = j;
        }
      }
    }
    parents[i] = smallest + dims[0] + 1;
  }
  for(int i = dims[0]; i < 2*dims[0] - 1; ++i) {
    smallest = dims[0] - 1;
    secsmallest = dims[0];
    for(int j = 0; j < dims[0] - 1; ++j) {
      cladeinclade = true;
      for(int k = 0; k < dims[0]; ++k) {
        if(bestcladetree[(i - dims[0]) * dims[0] + k]) {
          if(!bestcladetree[j * dims[0] + k]) {
            cladeinclade = false;
          }
        }
      }
      if(cladeinclade) {
        if(cladesizes[j] < cladesizes[smallest]) {
          secsmallest = smallest;
          smallest = j;
        } else if(cladesizes[j] < cladesizes[secsmallest]) {
          secsmallest = j;
        }
      }
    }
    if(secsmallest == dims[0] - 1) {
      parents[i] = 0;
    } else {
      parents[i] = secsmallest + dims[0] + 1;
    }
  }


  std::vector<int> result(4 * dims[0] - 2);
  for(int i = 0; i < dims[0]; ++i) {
    result[i] = parents[i];
    result[i + 2*dims[0] - 1] = dims[1];
  }
  for(int i = dims[0]; i < 2*dims[0] - 1; ++i) {
    result[i] = parents[i];
    result[i + 2*dims[0] - 1] = besttreescores[i - dims[0]];
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

