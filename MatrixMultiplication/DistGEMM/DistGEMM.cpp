#include "DistGEMM.hpp"
#include <MPI.h>
#include <libsci_acc.h>
#include <cassert>


DistGEMM::DistGEMM(int N, int P) {
	

	//build MPI topology
	int nums[3] = 
	
}

void DistGEMM::initializeLehmer() {

}

void DistGEMM::performGEMM() {
	// send root_k's matrix A to all other processors in communicator comm_k
	MPI_Bcast(A, blocksize*blocksize, mpi_val_type, root_k, comm_k);
	// send root_i's Matrix B to all other processors in communicator comm_i
	MPI_Bcast(B, blocksize*blocksize, mpi_val_type, root_i, comm_i);

	char transA = 'N';
	double alpha = 1.;
	double beta = 0.;

	dgemm(transA, transA, blocksize, blocksize, blocksize, alpha, A, N, B, N, beta, C, N);

	MPI_Reduce(C, C, blocksize*blocksize, mpi_val_type, MPI_SUM, root_j, comm_j);
}	
