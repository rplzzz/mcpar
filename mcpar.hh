#ifndef MCPAR_HH_
#define MCPAR_HH_

#include "mpi.h"
#include "mkl.h"
#include "vlfunc.hh"
#include "mcout.hh"
#include <stdlib.h>

class MCPar {
public:
  // Types, constants and enums
  enum {OK, INVALID, ERROR};
  
  const float TGT_ARATE_MIN;    // target acceptance rate minimum
  const float TGT_ARATE_MAX;    // target acceptance rate maximum
  const float SCALE_DEC;        // decrement factor when acceptance rate is too low
  const float SCALE_INC;        // increment factor when acceptance rate is too high
  const float PLOCAL;           // probablity of taking a local (instead of remote) proposal

  const int SYNCSTEP;           // frequency (number of steps) with
                                // which to synchronize gaussian
                                // posterior estimates

  static const float FPEPS;     // small number for preventing divide by zero and the like

  // logging switches.  Users may set and reset these as desired
  bool logging;
  int logstep;

  //! Constructor
  MCPar(int np, int nc=1, int mpisiz=1, int mpirank=0, float pl=0.9,
        float armin=0.2, float armax=0.5, float dfac=0.2, float ifac=1.5, int sync=10);
  ~MCPar();

  int run(int nsamp, int nburn, const float *pinit, VLFunc &L, MCout &outsamples,
          float *incov=0);
  void covar_setup(const float incov[], float *restrict cov);

  int genLocal(const float pvals[], float *restrict ptrial, float *restrict cfac);
  int genRemote(const float pvals[], float *restrict musigall, float *restrict ptrial,
                float *restrict cfac);
  // alternate genRemote for testing purposes only:
  int genRemote(const float pvals[], float * restrict ptrial, float * restrict cfac);

private:
  // communicator for use in the mcpar library
  MPI_Comm mcparComm;
  // parameters defining the problem size, number of parameters, etc. 
  int nparam;                   //<! number of parameters
  int nchain;                   //<! number of chains to run in parallel
  int ntot;                     //<! total number of parameters across
                                //! all the chains in this process
  int ncov;                     //! size of the covariance matrix (= nparam*nparam)

  // MPI configuration
  bool mpi;                     //<! flag indicating whether MPI is in use
  int rank;                     //<! MPI rank
  int size;                     //<! number of MPI processes
  int tchains;                  //<! total number of chains across all
                                //!processes

  // Working arrays for the MCMC

  // the pvals are ordered as {a1, b1, c1, a2, b2, c2, ...} for
  // parameters a,b,c and chains 1,2,... 
  float *pvals;                 // parameter values [ntot]
  float *ptrial;                // trial values for parameters [ntot]
  float *ptmp;                  // temporary storage for generating ptmp [ntot]
  float *mu;                    // parameter means, incrementally updated [ntot]
  float *sig;                   // parameter std. devs, incrementally updated [ntot]
  float *mutrial;               // temporary for mu in remote proposals [ntot]
  float *sigtrial;              // temporary for sigma^2 in remote proposals [ntot]
  float *qiarg;                 // temporary for use during the remote proposal evaluation [ntot]
  float *musigall;              // mu and sig from all chains across all processes [2*tchains*nparam]
  float *lylast;                // previous log-likelihood vals [nchain]
  float *lytrial;               // log-likelihood vals for ptrial [nchain]
  float *cfac;                  // Detailed balance correction factor for remote proposals [nchain]
  float *pacpt;                 // acceptance probability [nchain]
  float *acpt;                  // test for acceptance [nchain]
  float *qisum;                 // sum of the Q_i in remote proposals
  float *qimax;                 // max of the Q_i in remote proposals
  int *rjct;                    // continuation flag for rejection sampling [nchain]
  int *chnsel;                  // chain selector in remote proposals [nchain]

  // storage for the covariance matrix
  float *cov;                   // covariance matrix [ncov]

  // Arrays used in the running estimates of the local posterior.
  float *psum2;                 // numerator in variance calc.
  
  VSLStreamStatePtr rng;        // pointer to the rng 
};

#define VSL_CALL_CHK(c) {int stat = (c); if(stat != VSL_STATUS_OK) abort();}

#endif
