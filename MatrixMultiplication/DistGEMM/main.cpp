m#include <iostream>
#include <Measurement.hpp>
#include <Timer.hpp>
#include <cstdlib>
#include <omp.h>
#include <mpi.h>
#include <string>
#include <fstream>
#include <sstream>
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
	MPI_Barrier(MPI_COMM_WORLD);
//	Timer _t(1);
	double t1 = MPI_Wtime();
	m.performGEMM();
	MPI_Barrier(MPI_COMM_WORLD);
//	_t.stop();
	double t2 = MPI_Wtime();
	double t = t2-t1;
	Measurement mes("DistGEMM benchmark", N, N, t);

	int rank;
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank==0) {
		std::stringstream ss;
		ss << "benchmark-MPIWtime" << topsize << "x" << topsize << "x" << topsize << ".mes";	
		std::ofstream mesFile(ss.str().c_str(), std::ios::in | std::ios::app);
		mesFile << mes;
		mesFile.close();
	}
	//std::string filename("Msize=");
	//m.output_result(filename);
	MPI_Finalize();
	return 0;
}
