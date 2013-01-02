# Maurice Gonzenbach, 30.11.12, ETH Zurich
# compiles distributed gemm DistGEMM

CXX= CC
CXXFLAGS= -O3
INCLUDE=../HelperClasses/

LIBS=-lrt -openmp
#DEBUG= -g -DDEBUGOUTPUT -hlist=m 

all: DistGEMM

DistGEMM: main.o DistGemm.o
	$(CXX) $(DEBUG) $(CXXFLAGS) -I$(INCLUDE)  $(LIBS) main.o DistGemm.o -o DistGEMM.exe

main.o : main.cpp DistGEMM.cpp
	$(CXX) $(DEBUG) $(CXXFLAGS) -c -I$(INCLUDE) main.cpp -o main.o

DistGemm.o: DistGEMM.cpp DistGEMM.hpp
	$(CXX) $(DEBUG) $(CXXFLAGS) -c -I$(INCLUDE) $(LIBS) DistGEMM.cpp -o DistGemm.o

clean:
	rm -rf *.exe
	rm -rf *.o
	rm -r *.lst
