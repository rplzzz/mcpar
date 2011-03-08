#include "mkl.h"
#include "mkl_vsl.h"
#include "func.hh"
#include "mcpar.hh"

namespace mcp {

  enum {OK, INVALID, ERROR};
  
  int mcpar(int nsamp, int nburn, MCparm &parms, float *y, VLFunc &L, MCout &outsamples,
            MPIWrapper *mpi, float *incov)
  {
    int nparm  = parms.nparm();
    int nchain = parms.nchain();
    int ntot   = nparm*nchain;
    int rank   = mpi ? mpi.rank() : 0;
    int size   = mpi ? mpi.size() : 1; // number of MPI processes
    int tsize  = size*nchain;          // total number of chains across all processes
    

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
    float *pvals  = new float[ntot]; // parameter values
    float *ptrial = new float[ntot]; // trial parameters
    float *mu     = new float[ntot]; // parameter means - incremental update
    float *sig    = new float[ntot]; // parameter stdevs - incremental update
    floag *lylast = new float[nchain]; // previous log-likelihood values
    float *lytrial = new float[nchain]; // trial log-likelihood values
    float *pacpt  = new float[nchain]; // acceptance probability
    float *acpt   = new float[nchain]; // test for acceptance
    float scale   = 1.0f;
    float ntrial  = 0.0f;
    float naccept = 0.0f;
    
    // copy input guesses to our working array
    for(int j=0;j<ntot;++j)
      pvals[j] = parms.pvals[j];
    
    // first step - burn in.  Run the MC, but don't perform any remote updates
    int irate = nburn/10;
    for(int i = 0; i<nburn; ++i) {
      genLocal(nparm, nchain, parms, rng, ltrial);
      L(ltrial, ytrial);

      // XXX generate acpt from rng
      
      for(int j=0; j<nchain; ++j) {
        pacpt[j] = exp(ytrial[j]-ylast[j]);
        ylast[j] = acpt[j]<pacpt[j] ? ytrial[j] : ylast[j];
      }

      ntrial += nchain;
      for(int j=0; j<nchain; ++j) {
        int idx = j*nparm;
        if(acpt[j] < pacpt[j]) {
          naccept += 1.0f;
          for(int k=0; k<nparm; ++k,++idx)
            pvals[idx] = ptrial[idx];
        }
      }

      // tune acceptance rate
      if(i>irate) {
        float arate = naccept/ntrial;
        if(arate < TGT_ARATE_MIN) {
          naccept = ntrial = 0.0f; // reset the counters if we change size
          for(j=0;j<ncov;++j)
            cov[j] *= STEP_DEC;
        }
        else if(arate > TGT_ARATE_MAX) {
          naccept = ntrial = 0.0f; // reset counters
          for(j=0;j<ncov;++j)
            cov[j] *= STEP_INC;
        } 
        // evaluate again in 50 steps
        irate += 50;
      } 
    }

    // run the "keeper" samples, including remote proposals.  To keep
    // the calc well-vectorized, we do one check for whether to take a
    // local proposal or a remote one and apply it to all of the
    // chains in this process.
    for(int i=0; i<nsamp; ++i)
    
    // cleanup
    delete [] acpt;
    delete [] pacpt;
    delete [] sig;
    delete [] mu;
    delete [] ptrial;
    delete [] pvals;
    vsldeletestream(rng);
    delete [] cov;
  }
