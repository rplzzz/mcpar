#ifndef ROSENBROCK_HH_
#define ROSENBROCK_HH_

#include "vlfunc.hh"

/* Rosenbrock function with non-overlapping components */
class Rosenbrock1 : public VLFunc
{
private:
  const int n;                        // number of components.  Must be even

public:
  Rosenbrock1(int nc) : n(nc) {
    if(n<2 || n%2 != 0)
      throw("N for Rosenbrock1 must be even and >= 2"); // XXX replace with a real exception class
  }
  int operator()(int npset, const float *x, float *restrict fx);
};


/* Rosenbrock function with overlapping components */
class Rosenbrock2 : public VLFunc
{
private:
  const int n;                  // number of components. Must be >=2
public:
  Rosenbrock2(int nc) : n(nc) {
    if(n<2)
      throw("N for Rosenbrock2 must be >= 2"); // XXX replace with real exception class.
  }
  int operator()(int npset, const float *x, float *restrict fx);
};
    

/* Not really a Rosenbrock function at all; just a Gaussian */
class Gaussian : public VLFunc
{
private:
  const int n;                  // number of input variables -- fixed at 2 for testing
  float mu[2], sig2inv[2];
public:
  Gaussian(int nc, const float muin[]=0, const float sig2[]=0) : n(nc) {
    if(nc!=2) throw("Invalid specification.  N for Gaussian must == 2.");
    for(int i=0; i<n; ++i) {
      mu[i]      = muin ? muin[i] : 0.0f;
      sig2inv[i] = sig2 ? 1.0f/sig2[i] : 1.0f;
    } 
  }
  int operator()(int npset, const float *x, float *restrict fx);
};
  
#endif
