#ifndef MCOUT_HH_
#define MCOUT_HH_

#include <mpi.h>
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

  // members related to outputting results
  size_t nextout;               // offset of next parameter set to be output.  This is in terms of elements, not sets
  std::ostream *outstream;
  MPI_Comm mComm;
  int mrank;
  int msize;

public:
  MCout(int np, std::ostream *aoutstream, MPI_Comm acomm);
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
  void output();
  float * collect(size_t *ntot);
  void rewind(void) {nextout = 0;}
};


#endif
