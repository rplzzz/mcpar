
MKLBASE    = /cell_root/software/intel/mkl/10.0.2.018
MKLINC     = $(MKLROOT)/include
MKLLIBDIR  = $(MKLROOT)/lib/em64t
MKLLIBS    = -lmkl_solver_lp64_sequential -Wl,--start-group -lmkl_intel_lp64 -lmkl_sequential -lmkl_core -Wl,--end-group -lpthread
LIBS = #-lmpi

CXX      = mpiicpc
#CXX	  = mpicxx
#CXXFLAGS = -g -MMD -O0 -I$(MKLINC) -restrict
CXXFLAGS = -g -MMD -O3 -axSSE3 -fno-inline-functions -vec-report3 -I$(MKLINC) -restrict

OBJS = mcpar.o rosenbrock.o mcout.o mcutil.o
DEPS = $(OBJS:.o=.d)

include $(DEPS)

lib: libmcpar.a

all: mcpar-gauss mcpar-dgauss

libmcpar.a: $(OBJS)
	ar -ru libmcpar.a $(OBJS)

mcpar-gauss: mcpar-gauss.o $(OBJS)
	$(CXX) -L$(MKLLIBDIR) -o mcpar-gauss mcpar-gauss.o $(OBJS) $(LIBS) $(MKLLIBS)

mcpar-dgauss: mcpar-dgauss.o $(OBJS)
	$(CXX) -L$(MKLLIBDIR) -o mcpar-dgauss mcpar-dgauss.o $(OBJS) $(LIBS) $(MKLLIBS)

mcpar-dgauss-mpi: mcpar-dgauss-mpi.o $(OBJS)
	$(CXX) -L$(MKLLIBDIR) -o $@ $^ $(LIBS) $(MKLLIBS)

%.exe: %.o $(OBJS)
	$(CXX) -L$(MKLLIBDIR) -o $@ $^ $(LIBS) $(MKLLIBS)

clean:
	rm *.o
