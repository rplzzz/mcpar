#include "mcutil.hh"

void mcutil::qriguess(int rank, int npset, int nparam, const float plo[], const float phi[],
                      float * restrict pout)
{
  /* Initialize the qrng stream.  We do it here instead of in the
     constructor because we need to know nparam and don't want to
     provide it in the constructor.

     TODO: Fix this function so that it skips ahead by the number of
     MPI processes, so that we could call it multiple times and get
     new qrngs not used in any other process.
  */
  if(qrng)
    vslDeleteStream(&qrng); 
  VSL_CALL_CHK( vslNewStream(&qrng, VSL_BRNG_SOBOL, nparam) );

  
  int ntot = npset*nparam;
  float *qrngvals = new float[ntot];
  
  if(rank > 0)
    VSL_CALL_CHK( vslSkipAheadStream(qrng, rank*ntot) );
  
  VSL_CALL_CHK( vsRngUniform(VSL_RNG_METHOD_UNIFORM_STD, qrng, ntot, (float*)qrngvals,
                             0.0f, 1.0f) );

  for(int j=0; j<npset; ++j)
    for(int i=0; i<nparam; ++i) {
      int indx = j*nparam+i;
      pout[indx] = plo[i] + qrngvals[indx]*(phi[i]-plo[i]);
    }

}
