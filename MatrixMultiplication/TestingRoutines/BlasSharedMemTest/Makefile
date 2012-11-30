# Maurice Gonzenbach, 30.11.12, ETH Zurich
# compiles testing routines for different blas libraries
# USAGE: make PRGENV=compiler_used 

CXX= CC
CXXFLAGS= -std=c++11 -O3
INCLUDE=-I../../HelperClasses/
LIBS=

all: sci acml

sci: dgemm_libsci.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) dgemm_libsci.cpp $(LIBS) -lacml -o dgemm_libsci_$(PRGENV).exe

acml: dgemm_acml.cpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) dgemm_libsci.cpp $(LIBS) -lacml -o dgemm_libsci_$(PRGENV).exe

.PHONY:helpers
helpers:
	cd $(INCLUDE) && make

clean:
	rm *.exe && cd $(INCLUDE) && make clean
