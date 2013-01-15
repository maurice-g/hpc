#makefile to compile the 3ddiffusion code
#change next line to CC if on toedi, add -fast (-O3) to CXXFLAGS

CXX = mpicc
CXXFLAGS = #-O3 -hlist=m
INCLUDE = ./HelperClasses
LIBS = -lmpi -lmpi++

all: Diffusion3d

Diffusion3d: main.o Diffusion3d.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) main.o Diffusion3d.o -o Diffusion3d.exe

main.o : main.cpp Diffusion3D.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c  main.cpp -o main.o

Diffusion3d.o: Diffusion3D.cpp Diffusion3D.hpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c Diffusion3D.cpp -o Diffusion3d.o

clean:
	rm -rf *.exe
	rm -rf *.o

