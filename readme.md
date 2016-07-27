# MCPar:  A Parallel Markov Chain Monte Carlo Driver

## Introduction

MCPar is a parallel implementation of the Metropolis-Hastings Monte
Carlo algorithm (Metropolis, et al. 1953, Hastings 1970).  The
algorithm samples the distribution using multiple concurrent Markov
chains, correcting for the quasi-ergodicity problem, using a variant
of the algorithm proposed by Murray (2010).  (That paper also explains
the quasi-ergodicity problem, for readers unfamiliar with it.)  The
code is able to parallelize over multiple processors using MPI and
within a single processor using vectorization.

The code is designed so that users can plug in their own likelihood
functions.  Likelihood functions are implemented by writing subclasses
of the `VLFunc` interface class.  The version of the code in this
repository has a couple of simple examples using sums of Gaussians and
Rosenbrock functions (Rosenbrock, 1960).  There is also a wrapper for
a likelihood function written in R, if that is the sort of thing you
are into.

## Requirements

The code is set up to link to a set of MPI libraries.  It has been
tested with both OpenMPI and Intel MPI.  It also uses the Intel Math
Kernel Library (MKL).  You will need to have both of these libraries
in order to build and run the code.  

If you want to write your likelihood function in R, you will also need
to have R installed, along with the Rcpp and RInside packages.  This
functionality is optional and is not built by default.  If you want to
use it, set and export the environment variable USE_RFUNC before you
build.  

## References

Hastings, W. K. (1970), "Monte Carlo Sampling Methods Using Markov
Chains and Their Applications", _Biometrika_ **57**: 97--109

Metropolis, N., et al. (1953), "Equations of state calculations by
fast computing machines", _J. Chem. Phys._ **21**: 1087--92.

Murray, L. (2010), "Distributed Markov Chain Monte Carlo",
_Proceedings of Neural Information Processing Systems workshop on
learning on cores, clusters, and clouds_ **11**.

Rosenbrock, H. H. (1960), "An automatic method for finding the
greatest or least value of a function", _The Computer Journal_ **3**:
175--184.
