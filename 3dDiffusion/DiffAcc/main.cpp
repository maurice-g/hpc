//runfile for the 3dDiffusion simulation
//Stefan H & Maurice G
#include <iostream>
#include "Diffusion3D.hpp"
#include <mpi.h>
#include <string>
#include <omp.h>
#include <cstdlib>

int main(int argc, char *argv[]) {
	if (argc < 8) {
		std::cerr << "Too few arguments: [dx] [nx] [ny] [nz] [topology x] [topology y] [topology z]\n";
		return 1;
	}
	omp_set_num_threads(16);
	MPI_Init(&argc,&argv);
	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	
	double dx = atof(argv[1]);

	unsigned int nx = std::atoi(argv[2]);
	unsigned int ny = std::atoi(argv[3]);
	unsigned int nz = std::atoi(argv[4]);
	unsigned int topo_x = std::atoi(argv[5]);
	unsigned int topo_y = std::atoi(argv[6]);
	unsigned int topo_z = std::atoi(argv[7]);
	Diffusion3D::coord_type topology = {topo_x,topo_y,topo_z};	//a 2x2x2 mesh
	double D = 1.;
	double T = 0.4;


	Diffusion3D Simulation(dx,nx,ny,nz,topology,D,T);


	MPI_Comm_rank(MPI_COMM_WORLD,&rank);
	if (rank==0) {
		//Simulation.write_debug_info(std::cout);
	}
	//# flop: 9 per gridpoint
	MPI_Barrier(MPI_COMM_WORLD);
	double t0 = MPI_Wtime();
	Simulation.start_simulation(1);	
	MPI_Barrier(MPI_COMM_WORLD);
	double t1 = MPI_Wtime();
	double runtime = t1-t0;
	double flops = (nx*ny*nz*9*T/Simulation.get_dt())/runtime;
	std::string fname = "densities-";
//	Simulation.print_density(fname);
	if (rank==0) 
		std::cout << "3dDiffusion-nodes-dx-nx-ny-nz-runtime-flops " << topo_x*topo_y*topo_z << " " << dx << " " << nx << " " << ny << " " << nz <<
				  " " << runtime << " " << flops << "\n";
	MPI_Finalize();
	return 0;
}
