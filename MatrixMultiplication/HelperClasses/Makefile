# Maurice Gonzenbach, 30.11.12, ETH Zurich
# compiles the helper classes 

CXX= CC
CXXFLAGS=-O3
INCLUDE=
LIBS=-lrt

all: aligned_allocator.gch fillMatrix.gch matrix.gch Measurement.gch

%.gch: %.hpp
	$(CXX) $(CXXFLAGS) $(INCLUDE) $(LIBS) -c $<

clean:
	rm -rf *.gch 
