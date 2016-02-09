#include <Rcpp.h>
#include <vector>
#include <iterator>
#include <cmath>
#include "CCtrans_support.h"
using namespace Rcpp;

std::vector<int> ptroottrans(std::vector<int> pars, int ID) {
  std::vector<int> ans;
  ans.push_back(ID);
  int nextID = pars[ID - 1];
  while(nextID > 0) {
    ans.push_back(nextID);
    nextID = pars[nextID - 1];
  }
  return(ans);
}


void makecladearray(std::vector<bool> &ca, const std::vector<int> &pars, 
                    int n, int ss) {
  std::vector<int> parset(n);
  std::vector<int> pos;
  
  for(int i = 0; i < ss; ++i) {
    for(int ii = 0; ii< n; ++ii) {
      parset[ii] = pars[i*n + ii];
    }
    for(int j = 0; j < n; ++j) {
      pos = ptroottrans(parset, j + 1);
      for(unsigned int k = 0; k < pos.size(); ++k) {
        ca[i*n*n + (pos[k]-1)*n + j] = true;
      }
    }
  }
}


void tallyclades(std::vector<int> &pos, std::vector<int> &count,
                 std::vector<int> &refs, const std::vector<bool> &ca,
                 int n, int ss) {
  bool cladeequal = false;
  bool nodeequal = true;
  unsigned int nextclade = 0;
  int nextnode = 0;
  
  pos.push_back(0);
  count.push_back(1);
  
  for(int i = 1; i < n * ss; ++i) {
    while(!cladeequal && nextclade < pos.size()) {
      while(nodeequal && nextnode < n) {
        if(ca[i*n + nextnode] !=
           ca[pos[nextclade]*n + nextnode]) {
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
      pos.push_back(i);
      count.push_back(1);
      refs[i] = pos.size() - 1;
    } else {
      ++count[nextclade - 1];
      refs[i] = nextclade - 1;
    }
    cladeequal = false;
    nodeequal = true;
    nextclade = 0;
    nextnode = 0;
  }
}

void buildclades(std::vector<bool> &uc, std::vector<int> &cf, std::vector<int> &wc, 
                 int &nc, const std::vector<int> &pars, int n, int ss) {
  //make an array with all clades for all trees
  std::vector<bool> cladearray(n * n * ss);
  makecladearray(cladearray, pars, n, ss);
  
  //find positions of unique clades in the cladearray, plus clade frequency
  std::vector<int> uniquecladepos;
  tallyclades(uniquecladepos, cf, wc, cladearray, n, ss);
  
  //build the array with unique clades
  nc = cf.size();
  uc.resize(nc * n);
  for(int i = 0; i < nc; ++i) {
    for(int j = 0; j < n; ++j) {
      uc[n * i + j] = cladearray[n * uniquecladepos[i] + j];
    }
  }
}




void cladetimestats(std::vector<double> &sums, std::vector<double> &ssq,
                    const std::vector<double> &tims, 
                    const std::vector<int> &pc, 
                    const std::vector<int> &wc,
                    const std::vector<bool> &uc, int n, int ss) {

  for(int i = 0; i < n; ++i) { //each clade in best tree
    for(int j = 0; j < ss; ++j) { //each post tree
      for(int k = 0; k < n; ++k) { //each clade in post tree
        if(pc[i] == wc[j * n + k]) {          //if clades are equal
          sums[i] += tims[j * n + k];
          ssq[i] += tims[j * n + k] * tims[j * n + k];
        }
      }
    }
  }
}


//pc = postclade, will contain n positions of clades in 'uc' that are accepted 
//cr = claderoots, will contain n hosts which are the clade roots in the accepted clades (in that order)
//uc = array with all unique clades, with each n bool to indicate the hosts
//op = ordered pairs of clade indicators (positions in 'uc') and their frequencies (negative) among all posterior trees
//n = outbreak size, ss = sample size (number of trees)
void selectclades(std::vector<int> &pc, const std::vector<bool> &uc, 
                  const std::vector<std::pair<int, int> > &op, int n) {
  
  //during tree construction, freeelements will indicate for each accepted clade which hosts are not (yet) assigned
  //to subclades. In the end, each clade should still contain one of these: the clade root.
  std::vector<bool> freeelements(n * n);
  
  //the two most frequent clades are added to pc
  pc[0] = op[0].second;
  pc[1] = op[1].second;
  
  //for the first two clades, their freeelements are identified
  int clade0size = 0;
  int clade1size = 0;
  for(int i = 0; i < n; ++i) {
    clade0size += uc[pc[0] * n + i];
    clade1size += uc[pc[1] * n + i];
  }
  if(clade0size > clade1size) {
    for(int i = 0; i < n; ++i) {
      freeelements[i] = uc[pc[0] * n + i] && !uc[pc[1] * n + i];
      freeelements[i + n] = uc[pc[1] * n + i];
    }
  } else {
    for(int i = 0; i < n; ++i) {
      freeelements[i] = uc[pc[0] * n + i];
      freeelements[i + n] = !uc[pc[0] * n + i] && uc[pc[1] * n + i];
    }
  }
  
  
  //'nextclade' by 'nextclade' from 'op' is tested for compatibility with
  //accepted clades in 'pc'. If accepted, it is added to pc.
  int nextclade = 2;  //next clade to test
  int candclade = 0;  //the ID (in 'pc') of the candidate clade
  int cladecount = 2;  //nr of clades accepted thus far
  bool cladeaccept;  //whether to accept the clade undergoing testing
  bool overlap, candmore, accmore, anyfree;  //do clades overlap, 
  //does the candidate or accepted clade have non-overlapping hosts, 
  //will there remain free hosts in the accepted clade?
  std::vector<bool> candfree(n);  //which hosts will be free in the candidate clade?
  
  while(cladecount < n) {
    candclade = op[nextclade].second;
    for(int i = 0; i < n; ++i) {
      candfree[i] = uc[candclade * n + i];
    }
    cladeaccept = true;
    for(int i = 0; i < cladecount && cladeaccept; ++i) { //testing against all accepted clades
      overlap = false;
      candmore = false;
      accmore = false;
      anyfree = false;
      for(int j = 0; j < n; ++j) {
        if(uc[candclade * n + j] && uc[pc[i] * n + j]) {overlap = true;}
        if(uc[candclade * n + j] && !uc[pc[i] * n + j]) {candmore = true;}
        if(!uc[candclade * n + j] && uc[pc[i] * n + j]) {accmore = true;}
      }
      if(overlap && candmore && accmore) {  //candidate clade not nested -> not compatible
        cladeaccept = false;
      } else if(overlap && accmore) {   //accepted clade contains candidate clade
        for(int j = 0; j < n; ++j) {
          //will accepted clade keep at least one free host?
          if(!uc[candclade * n + j] && freeelements[i * n + j]) {anyfree = true;}
        }
        cladeaccept = anyfree;
      } else if(overlap && candmore) {   //candidate clade contains accepted clade
        for(int j = 0; j < n; ++j) {
          //which hosts in candidate clade will not be free?
          //and: will the candidate clade keep at least one free host?
          if(candfree[j]) {
            if(uc[pc[i] * n + j]) {
              candfree[j] = false;
            } else {anyfree = true;}
          }
        }
        cladeaccept = anyfree; 
      } else {cladeaccept = true;}   //candidate clade and accepted clade have different hosts
    }
    if(cladeaccept) {
      pc[cladecount] = candclade;  //add candidate to accepted clades
      for(int i = 0; i < cladecount; ++i) {  
        for(int j = 0; j < n; ++j) {
          freeelements[i * n + j] = freeelements[i * n + j] && !candfree[j]; //remove free hosts from older accepted clades
        }
      }
      for(int i = 0; i < n; ++i) {
        freeelements[cladecount * n + i] = candfree[i];  //register free hosts in new clade
      }
      ++cladecount;
    }
    ++nextclade;
  }
  
  //uc is ready, now determine root host in each clade
  std::vector<std::pair<int, int> > claderoots(n);
  for(int i = 0; i < n; ++i) {
    claderoots[i].second = pc[i];
    for(int j = 0; j < n; ++j) {
      if(freeelements[i * n + j]) {
        claderoots[i].first = j;
      }
    }
  }
  std::sort(claderoots.begin(), claderoots.end());
  
  for(int i = 0; i < n; ++i) {
    pc[i] = claderoots[i].second;
  }
  
}


//Function to find the parents given the bestcladetree parents, scores, bestcladetree, besttreescores, cladeparents
void findparents(std::vector<int> &pars, const std::vector<int> &pc, 
                 const std::vector<bool> &uc, int n) {
  std::vector<int> cladesizes(n + 1);
  for(int i = 0; i < n; ++i) {
    for(int j = 0; j < n; ++j) {
      cladesizes[i + 1] += uc[pc[i] * n + j];
    }
  }
  cladesizes[0] = n + 1;
  
  int parclade = 0;
  for(int i = 0; i < n; ++i) {
    parclade = 0;
    for(int j = 0; j < n; ++j) {
      if(uc[pc[j] * n + i] && !(j == i) && (cladesizes[j + 1] < cladesizes[parclade]) ) {
        parclade = j + 1;
      }
    }
    pars[i] = parclade;
  }
  
}

