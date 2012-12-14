# Maurice Gonzenbach, 30.11.12, ETH Zurich
# compiles distributed gemm DistGEMM

CXX= CC
CXXFLAGS= -O3
INCLUDE=../HelperClasses/

LIBS=-lrt -openmp

all: DistGEMM

DistGEMM: DistGEMM.hpp DistGEMM.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) main.cpp DistGEMM.cpp $(LIBS) -o DistGEMM.exe

clean:
	rm -rf *.exe
