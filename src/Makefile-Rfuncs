## -*- mode: make; tab-width: 8; -*-
##
## Simple Makefile
##
## TODO: 
##  proper configure for non-Debian file locations,   [ Done ]
##  allow RHOME to be set for non-default R etc

## comment this out if you need a different version of R, 
## and set set R_HOME accordingly as an environment variable
R_HOME := 		$(shell R RHOME)

sources := 		$(wildcard *.cpp)
programs := 		$(sources:.cpp=)


## include headers and libraries for R 
RCPPFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --cppflags)
RLDFLAGS := 		$(shell $(R_HOME)/bin/R CMD config --ldflags)
RBLAS := 		$(shell $(R_HOME)/bin/R CMD config BLAS_LIBS)
RLAPACK := 		$(shell $(R_HOME)/bin/R CMD config LAPACK_LIBS)

RRPATH :=		-Wl,-rpath,$(R_HOME)/lib

## include headers and libraries for Rcpp interface classes
## note that RCPPLIBS will be empty with Rcpp (>= 0.11.0) and can be omitted
RCPPINCL := 		$(shell echo 'Rcpp:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RCPPLIBS := 		$(shell echo 'Rcpp:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)


## include headers and libraries for RInside embedding classes
RINSIDEINCL := 		$(shell echo 'RInside:::CxxFlags()' | $(R_HOME)/bin/R --vanilla --slave)
RINSIDELIBS := 		$(shell echo 'RInside:::LdFlags()'  | $(R_HOME)/bin/R --vanilla --slave)

## compiler etc settings used in default make rules
CPPFLAGS += 		$(shell $(R_HOME)/bin/R CMD config CPPFLAGS)
CXXFLAGS += 		$(RCPPFLAGS) $(RCPPINCL) $(RINSIDEINCL) -DUSE_RFUNC
RFUNCLIBS := 		$(RLDFLAGS) $(RRPATH) $(RBLAS) $(RLAPACK) $(RCPPLIBS) $(RINSIDELIBS)

RFUNCOBJS = rfunc.o
RFUNCEXE  = mcpar-rfunc
RFUNCMAIN = mcpar-rfunc.o


all: $(RFUNCEXE)

rflgtest:
	@echo "RLDFLAGS    $(RLDFLAGS)"
	@echo "RRPATH      $(RRPATH)"
	@echo "RBLAS       $(RBLAS)"
	@echo "RLAPACK     $(RLAPACK)"
	@echo "RCPPLIBS    $(RCPPLIBS)"
	@echo "RINSIDELIBS $(RINSIDELIBS)"
	@echo "CXXFLAGS    $(CXXFLAGS)"

$(RFUNCEXE): $(OBJS) $(RFUNCOBJS) $(RFUNCMAIN)
	$(CXX) -L$(MKLLIBDIR) -o $@ $^ $(LIBS) $(RFUNCLIBS) $(MKLLIBS)
