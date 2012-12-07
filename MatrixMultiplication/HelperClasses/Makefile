# Maurice Gonzenbach, 30.11.12, ETH Zurich
# compiles the helper classes 

CXX= CC
CXXFLAGS= -std=c++11 -O3
INCLUDE=
LIBS=-lrt

all: aligned_allocator.o fillMatrix.o matrix.o Measurement.o

%.o: %.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $<

clean:
	rm -rf *.gch 
