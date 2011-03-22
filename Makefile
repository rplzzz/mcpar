
MKLBASE    = /cell_root/software/intel/mkl/10.0.2.018
MKLINC     = $(MKLROOT)/include
MKLLIBDIR  = $(MKLROOT)/lib/em64t
MKLLIBS    = -lmkl_solver_lp64_sequential -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
LIBS = #-lmpi

CXX      = icpc
CXXFLAGS = -g -I$(MKLINC) -restrict

OBJS = mcpar.o rosenbrock.o 

mcpar-gauss: mcpar-gauss.o $(OBJS)
	$(CXX) -L$(MKLLIBDIR) -o mcpar-gauss mcpar-gauss.o $(OBJS) $(LIBS) $(MKLLIBS)

