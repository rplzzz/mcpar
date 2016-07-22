#include <math.h>
#include "rosenbrock.hh"

int Rosenbrock1::operator()(int npset, const float *x, float *restrict fx)
{
  int ntot = npset * n;
  
  for(int j=0; j<npset; ++j)
    fx[j] = 0.0f;
  
  for(int i=0; i<ntot-1; i+=2) {
    int j=i/n;
    float t1 = 1-x[i];
    float t2 = x[i+1] - x[i]*x[i];
    
    // we define this with a negative sign so as to create a maximum
    // rather than a minimum.
    fx[j] -= t1*t1 + 100.0f*t2*t2;
  }
  return 0;
}



int Rosenbrock2::operator()(int npset, const float *x, float *restrict fx)
{
  int ntot = npset*n;

  for(int j=0; j<npset; ++j)
    fx[j] = 0.0f;

  for(int i=0; i<ntot-1; ++i) {
    int j=i/n;
    float t1 = 1-x[i];
    float t2 = x[i+1] - x[i]*x[i];

    // see note above
    fx[j] -= t1*t1 - 100.0f*t2*t2;
  }
  return 0;
}


int Gaussian::operator()(int npset, const float *x, float *restrict fx)
{
  int ntot = npset*n;

  for(int j=0; j<npset; ++j)
    fx[j] = 0.0f;

  for(int i=0; i<ntot; ++i) {
    int j=i/n;
    int k=i%n;
    float arg = x[i]-mu[k];
    fx[j] -= 0.5f*arg*arg*sig2inv[k];
  }

  // note that these goal functions are log-likelihood, so we don't
  // have to call exp here.
  return 0;
}

int DualGaussian::operator()(int npset, const float *x, float *restrict fx)
{
  const int n=2;                  // two parameters

  for(int j=0;j<npset;++j) {
    int indx = j*n;
    float arg1 = 0.5*(x[indx]*x[indx] + x[indx+1]*x[indx+1]);
    float t2a  = x[indx]-5.0f;
    float t2b  = x[indx+1]-5.0f;
    float arg2 = 0.5*(t2a*t2a + t2b*t2b);

    // have to return the log of sum of exponentials
    fx[j] = log(w*exp(-arg1) + exp(-arg2));
  }
  return 0;
}

  
