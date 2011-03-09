#include "mkl.h"
#include "mkl_vsl.h"
#include "func.hh"
#include "mcpar.hh"

  
int MCPar::run(int nsamp, int nburn, MCparm &parms, float *y, VLFunc &L, MCout &outsamples,
               MPIWrapper *mpi, float *incov)
{
  int nparm  = parms.nparm();   // XXX move to constructor
  int nchain = parms.nchain();  // XXX move to constructor
  int ntot   = nparm*nchain;
  int rank   = mpi ? mpi.rank() : 0;
  int size   = mpi ? mpi.size() : 1; // number of MPI processes
  int tsize  = size*nchain;          // total number of chains across all processes
  
  
  // set up the covariance matrix.  Copy incov if it was provided; otherwise use identity matrix
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
  
  
  // here are the things we will need to do the monte carlo
  float scale   = 1.0f;
  float ntrial  = 0.0f;
  float naccept = 0.0f;
    
  // copy input guesses to our working array
  for(int j=0;j<ntot;++j)
    pvals[j] = parms.pvals[j];
  
  // first step - burn in.  Run the MC, but don't perform any remote updates
  int irate = nburn/10;
  for(int i = 0; i<nburn; ++i) {
    genLocal(nparm, nchain, parms, rng, ptrial);
    L(ptrial, ytrial);
    
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
  for(int i=0; i<nsamp; ++i) {
    float rndlocl;            // random draw to see if we do local or remote proposal
    if(i<SYNCSTEP)
      rndlocl = 0.0f;        // no data for remote proposal yet.
    else
      //XXX get random from rng
      rndlocl = rndg();
    
    float cfac; // Metropolis-Hastings detailed balance correction factor
    if(rndlocl < PLOCAL) {
      genLocal(nparm, nchain, parms, rng, ptrial);
      lcfac = 1.0f;
    }
    else {
      genRemote(nparm, nchain, parms, muall, sigall, rng, ptrial, &cfac);
    }
    L(ptrial,ytrial);
    
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

  }
    
}

MCPar::MCPar(int np, int nc=1, int mpisiz=1, int mpirank=0) :
  nparm(np), nchain(nc), size(mpisiz), rank(mpirank)
{
  ntot = np*nc;
  ncov = np*np;

  tsize = mpisiz*nchain;
  mpi   = mpisiz>1;             // if mpisiz<=1, assume we're not running with MPI

  // allocate arrays
  pvals   = new float[ntot];
  ptrial  = new float[ntot];
  mu      = new float[ntot];
  sig     = new float[ntot];
  
  musig   = new float[2*tsize]; 

  lylast  = new float[nchain];
  lytrial = new float[nchain];
  pacpt   = new float[nchain];
  acpt    = new float[nchain];
  
  cov     = new float[ncov];

  // set up the rng
  // CAUTION!: the MKL only has 1024 independent MT rngs, so you can only do up to
  // 1024 processes at a time.  Within one process the random numbers are interleaved
  // amongst the markov chains, so you can have as many as you want.
  rng = NULL;
  if(rank>1024) {
    std::cerr << "Invalid rank > 1024 : rank = " << rank << "\n";
    throw("Invalid rank > 1024"); // XXX replace with a real exception class.
  } 
  int rngid = rank+VSL_BRNG_MT2203;
  vslNewStream(&rng, rngid, 8675309); 
}

MCPar::~MCPar()
{
  if(rng)
    vsldeletestream(rng);

  delete [] cov;
  delete [] acpt;
  delete [] pacpt;
  delete [] lytrial;
  delete [] lylast;
  delete [] musig;
  delete [] sig;
  delete [] mu;
  delete [] ptrial;
  delete [] pvals;
}
