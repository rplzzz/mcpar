# Example Job Submission Scripts

Since mcpar is designed to run on clusters, you will probably need to
run it by submitting a script to a batch queue.  This subdir contains
some examples that work with the Slurm batch scheduler.  _You can't
run these scripts directly._  They are meant to be submitted to a
batch scheduler using a command like `sbatch`.

By their nature these scripts are a little specific to the particular
problem and system, but these should give an idea of what you have to
do to run this code.  Some things users will definitely need to
customize:  

* The `#SBATCH -A` line specifies the account to charge the CPU time
  against.  You will need to replace this with the name of your
  account.

* The `module` commands load software packages into the environment.
  Your cluster may use different names for these modules, a different
  configuration manager, or no configuration manager at all.  Consult
  your system documentation to find out how to get the necessary
  software loaded in your environment.  You will need MPI, MKL, and a
  relatively recent compiler (from early 2016 or later).

* You will need to change file paths where applicable.

The job scripts provided here are:

* `mcpar-dgauss.sh`:  Run the Dual Gaussian demo calculation.

* `mcpar-rosen1.sh`:  Run the Rosenbrock function demo.

* `mcpar-rfunc.sh`: Run a calculation with a likelihood function
  written in R.  The particular R function and data are for another
  project and are not included in this repository, but you should be
  able to substitute your own likelihood function and input data set.

* `mcpar-fd.zsh`: Run an array of calculations with a different input
  data file for each calculation.  The file with the R code stays the
  same over all of the calculations in the array.  As before, the
  particular R model this script references isn't included here.  This
  script must be submitted with the `--array` option to `sbatch`.
