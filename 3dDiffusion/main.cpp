//runfile for the 3dDiffusion simulation
//Stefan H & Maurice G
#include <iostream>
#include "Diffusion3D.hpp"
#include <mpi.h>


int main(int argc, char ** argv) {
	MPI_Init(&argc,&argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0) {
		Vector3d<double> V(3,3,3);
		std::cout << V;
	}
	
	MPI_Finalize();
	return 0;
}
