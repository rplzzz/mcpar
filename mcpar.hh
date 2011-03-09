#ifndef MCPAR_HH_
#define MCPAR_HH_

class MCPar {
public:
  // Types and enums
  enum {OK, INVALID, ERROR};

  //! Constructor
  MCPar::MCPar(int np, int nc=1, int mpisiz=1, int mpirank=0);
  ~MCPar();
  
private:
  // parameters defining the problem size, number of parameters, etc. 
  int nparm;                    //<! number of parameters
  int nchain;                   //<! number of chains to run in parallel
  int ntot;                     //<! total number of parameters across
                                //! all the chains in this process
  int ncov;                     //! size of the covariance matrix (= nparm*nparm)

  // MPI configuration
  bool mpi;                     //<! flag indicating whether MPI is in use
  int rank;                     //<! MPI rank
  int size;                     //<! number of MPI processes
  int tsize;                    //<! total number of chains across all
                                //!processes

  // Working arrays.  We allocate these in the constructor, where
  // possible.  For the ones that depend on the MPI config, we
  // allocate them once we have the MPI info, but after that they are
  // not allowed to change size.  All the allocated parameters are
  // deallocated in the destructor.
  float *pvals;                 // parameter values [ntot]
  float *ptrial;                // trial values for parameters [ntot]
  float *mu;                    // parameter means, incrementally updated [ntot]
  float *sig;                   // parameter [ntot]
  float *musig;                 // mu and sig from all chains across all processes [2*tsize]
  float *lylast;                // previous log-likelihood vals [nchain]
  float *lytrial;               // log-likelihood vals for ptrial [nchain]
  float *pacpt;                 // acceptance probability [nchain]
  float *acpt;                  // test for acceptance [nchain]
  
  VSLStreamStatePtr rng;        // pointer to the rng
  

#endif
