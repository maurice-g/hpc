#include <iostream>
#include "DistGEMM.hpp"
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
int main(int argc, char* argv[]) {
	if (argc != 4) {
		std::cerr << "Usage: " << argv[0] << " [matrix dimension] [#procs] [3d topology size]" << std::endl;
		exit(-1);
	}
	MPI_Init(&argc,&argv);
	int N = std::atoi(argv[1]);
	int nprocs = std::atoi(argv[2]);
	int topsize = std::atoi(argv[3]);	//size of 3dtopology: topsize x topsize x topsize
	
	DistGEMM m(N,nprocs,topsize);
	m.initializeLehmer();
	m.performGEMM();
	//m.output_result();
	MPI_Finalize();
	return 0;
}
