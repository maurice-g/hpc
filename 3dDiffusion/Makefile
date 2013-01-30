#makefile to compile the 3ddiffusion code
#change next line to CC if on toedi, add -fast (-O3) to CXXFLAGS

CXX = CC
INTELFLAGS= -O3 -xHOST -ip -DNDEBUG -ansi-alias -fno-alias
CXXFLAGS =  -O3 -DNDEBUG -h restrict=a #-O3 -xHost -ip -DNDEBUG -ansi-alias -fno-alias #-hlist=m #-O3 -DNDEBUG -h restrict=a
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

D3AccCopy: main.o Diffusion3d_AccCopy.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) main.o Diffusion3d_AccCopy.o -o Diffusion3d.exe

Diffusion3d_AccCopy.o: Diffusion3D.hpp Diffusion3D_AccCopy.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c  Diffusion3D_AccCopy.cpp -o Diffusion3d_AccCopy.o
 
clean:
	rm -rf *.exe
	rm -rf *.o

