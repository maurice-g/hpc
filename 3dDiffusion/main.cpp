//runfile for the 3dDiffusion simulation
//Stefan H & Maurice G
#include <iostream>
#include "Diffusion3D.hpp"
#include <mpi.h>


int main(int argc, char ** argv) {
	MPI_Init(&argc,&argv);
	int rank, size;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	double dx = 0.1;
	unsigned int nx = 10;
	unsigned int ny = 10;
	unsigned int nz = 10;
	Diffusion3D::coord_type topology = {2,2,2};	//a 2x2x2 mesh
	double D = 1.;
	double T = 1.;

	
	Diffusion3D Simulation(dx,nx,ny,nz,topology,D,T);
	Simulation.write_debug_info(std::cout);
	
	
	MPI_Finalize();
	return 0;
}
