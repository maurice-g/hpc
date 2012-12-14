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
}	
