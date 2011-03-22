#ifndef MCOUT_HH_
#define MCOUT_HH_

#include <vector>

/* Class for storing output Monte Carlo samples.  Note the user is
 * responsible for ensuring that we don't overflow the available
 * storage.  Usually you will call newsamps(nsamp) right before
 * starting a loop that generates nsamp new samples. */
class MCout
{
  std::vector<float> pvals;
  const int nparam_;
  int next;
  int nset;                     // number of parameter sets stored

public:
  MCout(int np) : nparam_(np), next(0),nset(0) {}
  void newsamps(int nsamp) {
    int newsize = pvals.size() + nsamp*nparam_;
    pvals.resize(newsize);
  }
  void add(const float *pv) {
    float *strt = &pvals[next];
    for(int i=0; i<nparam_; ++i)
      strt[i] = pv[i];
    next += nparam_;
    nset++;
  }
  int size(void) {return nset;}
  int nparam(void) {return nparam_;}
  const float *pset(int i) const {return &pvals[i*nparam_];}
};


#endif
