#makefile to compile the 3ddiffusion code
#change next line to CC if on toedi, add -fast (-O3) to CXXFLAGS

CXX = CC
CXXFLAGS = -O3 -xHOST -ip -DNDEBUG -ansi-alias -fno-alias  #-hlist=m
INCLUDE = ./HelperClasses
LIBS = -openmp

all: Diffusion3d Diffusion3d-intel-optimized

Diffusion3d: main.o Diffusion3d.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) main.o Diffusion3d.o -o Diffusion3d.exe

Diffusion3d-intel-optimized: main.o Diffusion3d-opt.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) main.o Diffusion3d-opt.o -o Diffusion3d-opt.exe

main.o : main.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c  main.cpp -o main.o

Diffusion3d.o: Diffusion3D.cpp Diffusion3D.hpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c Diffusion3D.cpp -o Diffusion3d.o

Diffusion3d-opt.o: Diffusion3D-optimized.cpp Diffusion3D.hpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c Diffusion3D-optimized.cpp -o Diffusion3d-opt.o
clean:
	rm -rf *.exe
	rm -rf *.o

