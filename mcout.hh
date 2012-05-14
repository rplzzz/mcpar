#ifndef MCOUT_HH_
#define MCOUT_HH_

#include <vector>
#include <assert.h>
#include <iostream>

/* Class for storing output Monte Carlo samples.  Note the user is
 * responsible for ensuring that we don't overflow the available
 * storage.  Usually you will call newsamps(nsamp) right before
 * starting a loop that generates nsamp new samples. */
class MCout
{
  std::vector<float> pvals;
  const int nparam_;
  int next;
  int npset;                     // number of parameter sets stored
  int maxsamps;                  // maximum number of parameter sets that can be stored

public:
  MCout(int np) : nparam_(np), next(0),npset(0), maxsamps(0) {}
  void newsamps(int nsamp) {
    maxsamps += nsamp;
    int newsize = pvals.size() + nsamp*nparam_;
    pvals.resize(newsize);
  }
  void add(const float *pv) {
    float *strt = &pvals[next];
    for(int i=0; i<nparam_; ++i) {
      assert(strt+i < &pvals[0]+pvals.size());
      strt[i] = pv[i];
    }
    next += nparam_;
    npset++;
  }
  int size(void) const {return npset;}
  int maxsize(void) const {return maxsamps;}
  int nparam(void) {return nparam_;}
  int vsize(void) const {return pvals.size();}
  const float *getpset(int i) const {return &pvals[i*nparam_];}
  void output(std::ostream &outfile, int mpisize=1, int mpirank=0) const;
  float * collect(int mpisize, int mpirank, size_t *ntot) const;
};


#endif
