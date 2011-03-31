#include <iostream>
#include "mkl.h"
#include "mkl_vsl.h"
#include "mcpar.hh"
#include "vlfunc.hh"
#include "mcout.hh"
#include <math.h>

const float MCPar::FPEPS = 1.0e-14;

int MCPar::run(int nsamp, int nburn, const float *pinit, VLFunc &L, MCout &outsamples,
               float * incov)
{
  covar_setup(incov, cov);
  
  // pre-allocate space for the samples
  outsamples.newsamps(nsamp*nchain); 
  
  // here are the things we will need to do the monte carlo
  float ntrial  = 0.0f;
  float naccept = 0.0f;
  int stat;

  float *restrict pv = pvals;
  // copy input guesses to our working array
  for(int i=0;i<ntot;++i)
    pv[i] = pinit[i];

  // get the initial y values
  L(nchain, pvals, lylast);
  
  // first step - burn in.  Run the MC, but don't perform any remote updates
  int irate = 10;               // XXX reset to 50 when testing complete
  for(int isamp = 0; isamp<nburn; ++isamp) {
    genLocal(pvals, ptrial, cfac);
    L(nchain, ptrial, lytrial);

    // fill acpt with random numbers uniform on [0,1]
    VSL_CALL_CHK(vsRngUniform(VSL_METHOD_SUNIFORM_STD, rng, nchain, acpt, 0.0f, 1.0f)) ;
    
    ntrial += nchain;
    for(int j=0; j<nchain; ++j) {
      pacpt[j] = exp(lytrial[j]-lylast[j]);
      lylast[j] = acpt[j]<pacpt[j] ? lytrial[j] : lylast[j];
      naccept  += acpt[j]<pacpt[j];
    } 
    for(int i=0; i<ntot; ++i) {
      int j = i/nparam;
      if(acpt[j] < pacpt[j])
        pvals[i] = ptrial[i];
    }

    // tune acceptance rate
    if(isamp>irate) {
      float arate = naccept/ntrial;
      if(arate < TGT_ARATE_MIN) {
        naccept = ntrial = 0.0f; // reset the counters if we change size
        // scale covariance matrix.
        // XXX cov actually contains the Cholesky factorization of the covariance
        //     matrix.  Should we therefore be multiplying by the sqrt of the factor
        //     we really want
        for(int i=0; i<ncov; ++i)
          cov[i] *= SCALE_DEC;  // XXX some wasted effort here, since cov is symmetric
      }
      else if(arate > TGT_ARATE_MAX) {
        naccept = ntrial = 0.0f; // reset counters
        for(int i=0; i<ncov; ++i)
          cov[i] *= SCALE_INC;  // XXX wasted effort.
      } 
      // evaluate again in 50 steps
      irate += 10;              // XXX reset to 50 when testing complete
    } 
  }

  // initialize mu and psum2
  for(int i=0;i<ntot;++i) {
    mu[i] = 0.0f;
    psum2[i] = FPEPS;           // prevent divide by zero
  }
  float pwgt=0.0f;               // total weight in the mu, sigma calcs.
  
  // run the "keeper" samples, including remote proposals.  To keep
  // the calc well-vectorized, we do one check for whether to take a
  // local proposal or a remote one and apply it to all of the
  // chains in this process.
  for(int isamp=0; isamp<nsamp; ++isamp) {
    float rndlocal;            // random draw to see if we do local or remote proposal
    if(isamp<SYNCSTEP)
      rndlocal = 0.0f;        // no data for remote proposal yet.
    else
      VSL_CALL_CHK(vsRngUniform(VSL_METHOD_SUNIFORM_STD, rng, 1, &rndlocal, 0.0f, 1.0f)) ;
    
    // Should cfac include the total pdf (i.e., cfac = 0.9*cflocal +
    // 0.1*cfremote, irrespective of whether we are doing a remote or
    // local update)?
    int remotep;                // was this a remote proposal?
    if(rndlocal <= PLOCAL) {     
      genLocal(pvals, ptrial, cfac);
      remotep = 0;
    }
    else {
      genRemote(pvals, musigall, ptrial, cfac);
      //genRemote(pvals, ptrial, cfac);
      remotep = 1;
    }
    L(nchain, ptrial,lytrial);
    
    // fill acpt with random numbers uniform on [0,1]
    VSL_CALL_CHK(vsRngUniform(VSL_METHOD_SUNIFORM_STD, rng, nchain, acpt, 0.0f, 1.0f)) ;
    
    ntrial += nchain;
    for(int j=0; j<nchain; ++j) {
      pacpt[j] = exp(lytrial[j]-lylast[j]) * cfac[j];
      lylast[j] = acpt[j]<pacpt[j] ? lytrial[j] : lylast[j];
      naccept  += acpt[j]<pacpt[j];
    } 
    for(int i=0; i<ntot; ++i) {
      int j = i/nparam;
      if(acpt[j] < pacpt[j])
        pvals[i] = ptrial[i];
    }
    // Add the new samples to the list
    for(int j=0; j<nchain; ++j) {
      int idx = j*nparam;
      outsamples.add(pvals+idx);
    } 

    // Update the running estimates of mean and variance for each
    // chain
    pwgt += 1.0f;
    float winv = 1.0f/pwgt;
    for(int i=0; i<ntot; ++i) {
      int j     = i/nparam;
      if(remotep && acpt[j] < pacpt[j]) { // accepted a remote proposal
        // replace mu and sigma with the remote values.  We have to
        // reconstruct psum2.  Since we're plugging in (effectively)
        // *last* iteration's sigma^2, we do that by multiplying by n-1.
        mu[i]    = mutrial[i];
        sig[i]   = sigtrial[i];
        psum2[i] = sig[i]*(pwgt-1.0f);
      }
      
      float delta = pvals[i] - mu[i];
      mu[i] += delta * winv;
      psum2[i] += delta*(pvals[i]-mu[i]); // one factor of old mean, one factor of new
      sig[i] = psum2[i] * winv; // technically should be /(n-1), but
                                // it's not worth the extra fdiv.
      
      // Write the new estimates of mu and sigma into the MPI array.
      int islot = 2*(rank*ntot+i);
      musigall[islot]   = mu[i];
      musigall[islot+1] = sig[i];
    }                                    /* end of loop over all params */

    // sync with MPI peers, if any
    if(mpi && (isamp % SYNCSTEP == 0)) {
      // do MPI_Allgather:  not yet implemented 
    } 
  }                                      /* end of loop over samples */
  return 0;
}

MCPar::MCPar(int np, int nc, int mpisiz, int mpirank, float pl,
             float armin, float armax, float dfac, float ifac, int sync) :
  nparam(np), nchain(nc), size(mpisiz), rank(mpirank), PLOCAL(pl),
  TGT_ARATE_MIN(armin), TGT_ARATE_MAX(armax),
  SCALE_DEC(dfac), SCALE_INC(ifac), SYNCSTEP(sync)
{
  ntot = np*nc;
  ncov = np*np;

  tchains = mpisiz*nchain;
  mpi   = mpisiz>1;             // if mpisiz<=1, assume we're not running with MPI

  // allocate arrays
  pvals   = new float[ntot];
  ptrial  = new float[ntot];
  ptmp    = new float[ntot];
  mu      = new float[ntot];
  sig     = new float[ntot];
  mutrial = new float[ntot];
  sigtrial= new float[ntot];
  qiarg   = new float[ntot];
  
  musigall= new float[2*tchains*nparam]; 

  lylast  = new float[nchain];
  lytrial = new float[nchain];
  cfac    = new float[nchain];
  pacpt   = new float[nchain];
  acpt    = new float[nchain];
  qisum   = new float[nchain];
  qimax   = new float[nchain];
  rjct    = new int[nchain];
  chnsel  = new int[nchain];
  
  cov     = new float[ncov];

  psum2   = new float[ntot];

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
    vslDeleteStream(&rng);

  delete [] psum2;
  delete [] cov;
  delete [] chnsel;
  delete [] rjct;
  delete [] qimax;
  delete [] qisum;
  delete [] acpt;
  delete [] pacpt;
  delete [] cfac;
  delete [] lytrial;
  delete [] lylast;
  delete [] musigall;
  delete [] qiarg;
  delete [] sigtrial;
  delete [] mutrial;
  delete [] sig;
  delete [] mu;
  delete [] ptmp;
  delete [] ptrial;
  delete [] pvals;
}


int MCPar::genLocal(const float pvals[], float *restrict ptrial, float *restrict cfac)
{
  for(int j=0; j<nchain; ++j) {
    int indx = j*nparam;
    VSL_CALL_CHK(vsRngGaussianMV(VSL_METHOD_SGAUSSIANMV_BOXMULLER2, rng, 1,
                                 ptrial+indx, nparam, VSL_MATRIX_STORAGE_FULL,
                                 pvals+indx, cov ));
    cfac[j] = 1.0f;             // gaussian is symmetric
  }
  return 0;
}


int MCPar::genRemote(const float pvals[], float * restrict musigall,
                     float * restrict ptrial, float * restrict cfac)
{

  const float mutest[] = {1.0f,2.0f, 0.0f,2.0f, 0.0f,2.0f, 1.0f,2.0f, -1.0f,2.0f,
                          0.0f,2.0f, 0.0f,2.0f, -1.0f,2.0f};
  
//   for(int i=0;i<16;++i)
//     musigall[i] = mutest[i];

  
  // Each mu, sigma set (one pair for each parameter) defines a
  // distribution Q_i.  We want to choose over the distribution
  // max_i(Q_i), and we will do that by choosing from Sum_i(Q_i) and
  // doing rejection sampling.  Note that choosing from Sum_i(Q_i) is
  // the same as selecting an i with uniform probability and choosing
  // from Q_i.

  for(int j=0;j<nchain;++j)
    rjct[j] = 1;            // rjct==1 means continue sampling
  int anyrjct = 1;

  do {
  
    // get random integer from 0..tchains-1 for each chain.  We will use
    // that to select a Q_i.
    int stat = viRngUniform(VSL_METHOD_SUNIFORM_STD, rng, nchain, chnsel, 0, tchains);
    
    for(int j=0; j<nchain; ++j) {
      if(rjct[j]) {
        for(int i=0; i<nparam; ++i) {
          int tindx = 2*(nparam*chnsel[j] + i); // "total index" 
          int cindx  = j*nparam+i;                // "chain index"
          mutrial[cindx]  = musigall[tindx]; 
          // recall "sig" is actually sig^2 -- have to take sqrt to get Cholesky decomp
          sigtrial[cindx] = sqrt(musigall[tindx+1]);
        } 
        VSL_CALL_CHK(vsRngGaussianMV(VSL_METHOD_SGAUSSIANMV_BOXMULLER2, rng, 1,
                                     ptrial+j*nparam, nparam, VSL_MATRIX_STORAGE_DIAGONAL,
                                     mutrial+j*nparam, sigtrial+j*nparam )) ;
      }
    }
    

    for(int j=0; j<nchain; ++j)
      if(rjct[j]) {
        qimax[j] = FPEPS;
        qisum[j] = FPEPS;
      }
      else {                    // new params have been accepted for this chain
        // ensure that we have a zero pacpt at the end, so we don't
        // disturb the values accepted by this chain.
        qimax[j] = 0.0f;
        qisum[j] = 1.0f;
      }
    
    for(int qi=0; qi<tchains; ++qi) { // loop over all Q_i
      for(int i=0; i<ntot; ++i) {     // all parameters in all chains
        int indx  = 2*(qi*nparam + i%nparam); // index in the master mu-sigma array
        float xm   = musigall[indx] - ptrial[i];  // x-mu
        float sig2 = musigall[indx+1];            // sig^2
        
        qiarg[i] = xm*xm/sig2;
      }
      
      // qiarg now has the (x-mu)/sig2 terms for this Q_i.  
      for(int j=0; j<nchain; ++j) { // loop over chains
        float arg = 0.0f;

        // sum up all the qiargs corresponding to this chain
        for(int i=0; i<nparam; ++i) {
          int indx = j*nparam+i; // indx runs from 0 to ntot-1 on successive iterations of the j loop
          arg += qiarg[indx];
        }

        float gv = exp(-0.5*arg);   // Gaussian value

        // skip the final assignment, if this chain has already accepted a new
        // set of parameters.
        if(rjct[j]) {
          qisum[j] += gv;
          qimax[j]  = gv > qimax[j] ? gv : qimax[j];
        }
      }
    }
    
    for(int j=0; j<nchain; ++j)
      pacpt[j] = qimax[j]/qisum[j];      // ==0.0 for any chains that have already been accepted
    
    // generate acpt
    VSL_CALL_CHK(vsRngUniform(VSL_METHOD_SUNIFORM_STD, rng, nchain, acpt, 0.0f, 1.0f));

    anyrjct = 0;
    
    for(int j=0; j<nchain; ++j) {
      //if(acpt[j] < pacpt[j]) {  // accept this chain
      if(1) {
        rjct[j] = 0;            // mark as accepted -- won't get any
                                // further changes even if we go
                                // through the loop again
        
        
        // Formula for cfac is: cfac[j] = max_i(Q_i(pvals)) / max_i(Q_i(ptrial))
        // qimax[j] has the max for the ptrial values, so we just need the pvals
        // version
        
        // NB: The structural difference between this loop and the one
        // for ptrial is unsatisfying.  We should decide which way of
        // doing this calculation is better and stick with it.  We may
        // also want to move it to a subroutine so that it can be used
        // in both places.
        cfac[j]  = 0.0f;
        int ip0  = j*nparam;
        int ip1  = (j+1)*nparam;
        for(int qi=0; qi<tchains; ++qi) { // loop over all Q_i
          float arg= 0.0f;
          for(int i=ip0; i<ip1; ++i) { // all the parameters in this chain.
            int indx  = 2*(qi*nparam + i%nparam);
            float xm   = musigall[indx] - pvals[i];   // x-mu
            float sig2 = musigall[indx+1];            // sig^2
            
            arg += xm*xm/sig2;
          }

          // arg is the sum over all Q_i of (x-mu_i)^2/sig^2
          float gv = exp(-0.5*arg); // Gaussian value for Q_i   // XXX Move me outside the loop!
          //cfac[j] = gv > cfac[j] ? gv : cfac[j];
          cfac[j] += gv;        // XXX fix when testing complete
        }

        //cfac[j] /= qimax[j];  // XXX fix when testing complete.
        cfac[j] /= qisum[j];
      } /* end of if acpt[j]<pacpt[j] */
      anyrjct += rjct[j];       // do this for all j -- will be nonzero if any samples were rejected.
    } 
  } while(anyrjct);             // finish only when no rejects left.

  // sigtrial was set to 1/sigma^2 for use in the gaussian rng.
  // Invert it back so that we can resume our incremental updates.
//   for(int i=0;i<ntot;++i)
//     sigtrial[i] = 1.0f/sigtrial[i];
  
  return 0;                     // haven't created any return codes
}


void MCPar::covar_setup(const float *incov, float *restrict cov)
{
  // set up the covariance matrix.  Copy incov if it was provided; otherwise use identity matrix
  if(incov)
    for(int i=0; i<ncov; ++i)
      cov[i] = incov[i];
  else {
    for(int i=0; i<ncov; ++i)
      cov[i] = 0.0f;
    for(int i=0; i<nparam; ++i) {
      int indx = i*(nparam+1);  // ith diagonal element
      cov[indx] = 1.0f;
    }
  }

  // vsRngGaussianMV wants the Cholesky factorization of the
  // covariance matrix.  Much of what follows is boilerplate from the
  // Intel MKL examples

  
  /* Variables needed for Cholesky factorization */
  char uplo='U';
  int n = nparam;
  int lda = nparam;
  int stat;
  
  spotrf(&uplo, &n, cov, &lda, &stat);

  // Modulo some scaling, this factorization should be good for the
  // remainder of the run.
}

  
int MCPar::genRemote(const float pvals[], float * restrict ptrial, float * restrict cfac)
{
  // pick from one of 4 fixed gaussian distributions
  const float muvals[8] = {1.0f,0.0f, 0.0f,1.0f, -1.0f,0.0f, 0.0f,-1.0f};
  // const float muvals[8] = {0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f, 0.0f,0.0f};
  float unitcov[2] = {2.0f,2.0f};
  float uc2[2];
  float fullcov[4] = {0.5f, 0.0f, 0.0f, 0.5f}; 
  float ucmat[4];
  const int ndist = 4;
  const int gaussian = 1;

  for(int i=0;i<2;++i)
    uc2[i] = unitcov[i]*unitcov[i];
  
  covar_setup(fullcov,ucmat);
  
  int stat = viRngUniform(VSL_METHOD_SUNIFORM_STD, rng, nchain, chnsel, 0, ndist);

  for(int j=0;j<nchain;++j) {
    int muindx = 2*chnsel[j];

    if(gaussian) {
      VSL_CALL_CHK(vsRngGaussianMV(VSL_METHOD_SGAUSSIANMV_BOXMULLER2, rng, 1,
                                   ptrial+j*nparam, nparam, VSL_MATRIX_STORAGE_DIAGONAL,
                                   muvals+muindx, unitcov) );
//       VSL_CALL_CHK(vsRngGaussianMV(VSL_METHOD_SGAUSSIANMV_BOXMULLER2, rng, 1,
//                                    ptrial+j*nparam, nparam, VSL_MATRIX_STORAGE_FULL,
//                                    muvals+muindx, ucmat) );
    }
    else {
      VSL_CALL_CHK(vsRngUniform(VSL_METHOD_SUNIFORM_STD, rng, nparam, ptrial+j*nparam, -5.0f, 5.0f));
    }
  }
  
  // calculate cfac
  // cfac is Q(pvals)/Q(ptrial)
  for(int j=0; j<nchain; ++j) {
    int pindx  = nparam*j;
    float cpv  = 0.0f;
    float cpt  = 0.0f;
    
    for(int dist=0; dist<ndist; ++dist) {
      int muindx = 2*dist;
      float arg1 = pvals[pindx]   - muvals[muindx];
      float arg2 = pvals[pindx+1] - muvals[muindx+1];
      cpv += exp(-0.5*(arg1*arg1/uc2[0] + arg2*arg2/uc2[1]));

      arg1 = ptrial[pindx]   - muvals[muindx];
      arg2 = ptrial[pindx+1] - muvals[muindx+1];
      cpt += exp(-0.5*(arg1*arg1/uc2[0] + arg2*arg2/uc2[1])); 
    }

    if(gaussian)
      cfac[j] = cpv/cpt;
    else
      cfac[j] = 1.0f;
  }
  return 0;
}
