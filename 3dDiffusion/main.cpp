//runfile for the 3dDiffusion simulation
//Stefan H & Maurice G
#include <iostream>
#include "Diffusion3D.hpp"
#include <mpi.h>


int main(int argc, char ** argv) {
	MPI_Init(&argc,&argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double dx = 0.25;
	unsigned int nx = 2;
	unsigned int ny = 4;
	unsigned int nz = 6;
	Diffusion3D::coord_type topology = {2,2,2};	//a 2x2x2 mesh
	double D = 1.;
	double T = 1.;

	
	Diffusion3D Simulation(dx,nx,ny,nz,topology,D,T);
	

	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0) {
		Simulation.write_debug_info(std::cout);
	}
	
	
	MPI_Finalize();
	return 0;
}
