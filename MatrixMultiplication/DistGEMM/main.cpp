#include <iostream>
#include <Measurement.hpp>
#include <Timer.hpp>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
#include <string>
#include "DistGEMM.hpp"

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

	std::string matrixA("matrixA.in");
	std::string matrixB("matrixB.in");
	
	//m.setup(matrixA, matrixB);
	m.initializeLehmer();
	Timer _t(1);
	m.performGEMM();
	MPI_Barrier();
	_t.stop();
	Measurement mes("Distribution & DGEMM (libsci) operation", N, N, _t.elapsed_s());
	std::cout << mes; 
	//std::string filename("Msize=");
	//m.output_result(filename);
	MPI_Finalize();
	return 0;
}
