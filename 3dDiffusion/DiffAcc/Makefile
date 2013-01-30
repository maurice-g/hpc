#makefile to compile the 3ddiffusion code
#change next line to CC if on toedi, add -fast (-O3) to CXXFLAGS

CXX = CC
INTELFLAGS= -O3 -xHOST -ip -DNDEBUG -ansi-alias -fno-alias
CXXFLAGS =   -hlist=m  -O3 -DNDEBUG #-h restrict=a
INCLUDE = ./HelperClasses
LIBS = -openmp

all: Diffusion3d Diffusion3dCuda

Diffusion3d: main.o Diffusion3d.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) main.o Diffusion3d.o -o Diffusion3d.exe

Diffusion3dCuda: main.o Diffusion3dCuda.o cudaFunctions.o
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) main.o cudaFunctions.o Diffusion3dCuda.o -o Diffusion3dCuda.exe

main.o : main.cpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c  main.cpp -o main.o

Diffusion3d.o: Diffusion3D.cpp Diffusion3D.hpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c Diffusion3D.cpp -o Diffusion3d.o

Diffusion3dCuda.o: Diffusion3DCuda.cpp Diffusion3D.hpp
	$(CXX) $(CXXFLAGS) -I$(INCLUDE) $(LIBS) -c Diffusion3DCuda.cpp -o Diffusion3dCuda.o

cudaFunctions.o: cudaFunctions.cu
	nvcc  -c cudaFunctions.cu -o cudaFunctions.o
 
clean:
	rm -rf *.exe
	rm -rf *.o
	rm -rf *.lst
	rm -rf *.ptx
	rm -rf *.cub
