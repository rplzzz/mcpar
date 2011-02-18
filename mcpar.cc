#include "mkl.h"
#include "mkl_vsl.h"
#include "func.hh"
#include "mcpar.hh"

namespace mcp {

  enum {OK, INVALID, ERROR};
  
  int mcpar(MCparm &parms, float *y, VLFunc &L, MCout &outsamples,
            MPIWrapper *mpi, float *incov)
  {
    int nparm  = parms.nparm();
    int nchain = parms.nchain();
    int rank   = mpi ? mpi.rank() : 0;

    // set up the covariance matrix.  Copy incov if it was provided; otherwise use identity matrix
    int ncov   = nparam*nparam;
    float *cov = new float[ncov];
    if(incov)
      for(int i=0; i<ncov; ++i)
        cov[i] = incov[i];
    else {
      for(int i=0; i<ncov; ++i)
        cov[i] = 0.0f;
      for(int i=0; i<nparam; ++i) {
        int indx = i*(nparam+1);
        cov[indx] = 1.0f;
      }
    }

    // set up the rng
    // CAUTION!: the MKL only has 1024 independent MT rngs, so you can only do up to
    // 1024 processes at a time.  Within one process the random numbers are interleaved
    // amongst the markov chains, so you can have as many as you want.
    if(rank>1024) {
      std::cerr << "Invalid rank > 1024 : rank = " << rank << "\n";
      return ERROR;
    }
    int rngid = rank+VSL_BRNG_MT2203;
    VSLStreamStatePtr rng;
    vslNewStream(&rng, rngid, 8675309);
        

    // here are the things we will need to do the monte carlo
    float *trial = new float[nparm]; // trial parameters
    float *mu    = new float[nparm]; // parameter means - incremental update
    float *sig   = new float[nparm]; // parameter stdevs - incremental update
    float scale  = 1.0f;
    int   ntrial = 0;
    int  naccept = 0;
    

    
    
    // cleanup
    vsldeletestream(rng);
    delete [] cov;
  }
