#ifndef MCUTIL_HH_
#define MCUTIL_HH_

#include <mkl.h>
#include <mkl_vsl.h>


#define VSL_CALL_CHK(c) {int stat = (c); if(stat != VSL_STATUS_OK) abort();}


class mcutil {                  // various utility functions for mcmc

  VSLStreamStatePtr qrng;       // pointer to the qrng

public:
  mcutil() : qrng(0) {}
  ~mcutil() {if(qrng) vslDeleteStream(&qrng);}

  
  /*!
   * \brief quasi-random initial guess
   * \param[in] rank MPI rank of this process (=0 for serial mcmc)
   * \param[in] npset Number of parameter sets in this process
   * \param[in] nparam Number of parameters
   * \param[in] plo Lower bound for initial parameter guesses (float plo[nparam])
   * \param[in] phi Upper bound for initial parameter guesses (float phi[nparam])
   * \param[out] pout Output parameter sets (float pout[nparam*npset])
   */
  void qriguess(int rank, int npset, int nparam, const float plo[], const float phi[], float * restrict pout);

};

#endif
