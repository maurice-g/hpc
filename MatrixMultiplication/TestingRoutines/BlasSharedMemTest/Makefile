# Maurice Gonzenbach, 30.11.12, ETH Zurich
# compiles testing routines for different blas libraries
# USAGE: make PRGENV=compiler_used 

CXX= CC
CXXFLAGS= -O3
INCLUDE=../../HelperClasses/

ifeq ($(PRGENV),gnu)
	LOPENMP=-fopenmp
else
	LOPENMP=-openmp
endif

LIBS=$(LOPENMP) -lrt

all: libsci acml magma cublas cublas_pinned

libsci: dgemm_libsci.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) dgemm_libsci.cpp $(LIBS) -o dgemm_libsci_$(PRGENV).exe

acml: dgemm_acml.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) dgemm_acml.cpp $(LIBS) -o dgemm_acml_$(PRGENV).exe

magma: dgemm_magma.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) -I../../../../magma-1.3.0/include/ -L../../../../magma-1.3.0/lib/ dgemm_magma.cpp -L$(LIBS) -lmagma -lmagmablas -lcublas -fopenmp -o dgemm_magma_$(PRGENV).exe

cublas: dgemm_cublas.cu
	nvcc $(CXXFLAGS) -I$(INCLUDE) -lcublas dgemm_cublas.cu -o dgemm_cublas_$(PRGENV).exe -lrt

cublas_pinned: dgemm_cublas_hostalloc.cu
	nvcc $(CXXFLAGS) -I$(INCLUDE) -lcublas dgemm_cublas_hostalloc.cu -o dgemm_cublas_hostalloc_$(PRGENV).exe -lrt

#.PHONY:helpers
#helpers:
#	cd $(INCLUDE) && make

clean: clean-headers
	rm -rf *.exe

clean-headers:
	cd $(INCLUDE) && make clean
