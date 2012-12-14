#include "DistGEMM.hpp"
#include <MPI.h>
#include <libsci_acc.h>
#include <cassert>
#include <cmath>

DistGEMM::DistGEMM(int n, int numprocs, int cubes) {
	std::assert(cubes*cubes*cubes == numprocs);	
	std::assert(n % cubes == 0);
	
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	MPI_Comm_rank(MPI_COMM_WORLD,&rank);

	P = numprocs;
	N = n;
	cubeSize = cubes;					//# processors in each direction of cartesian topology

	//build MPI topology
	int nums[3];						//number of nodes in each dimension (x,y,z)
	MPI_Dims_create(size,3,nums);
	
	//build global communicator, then subdivide into comm_i, comm_j, comm_k
	int periodic[3] = {0, 0, 0};
	MPI_Comm cart_comm;
	MPI_Cart_create(MPI_COMM_WORLD,3,nums,periodic,true,&cart_comm);
	
	//assign comm_i, comm_j, comm_k
	int logv_i = {true,false,false};
	int logv_j = {false,true,false};
	int logv_k = {false,false,true};
	MPI_Cart_sub(cart_comm,logv_i,comm_i);
	MPI_Cart_sub(cart_comm,logv_j,comm_j);
	MPI_Cart_sub(cart_comm,logv_k,comm_k);
	
	//assign the coordinates to each rank
	int d = 3;
	int coords[3];
	MPI_Cart_coords(cart_comm,rank,d,coords);
	p_i = coords[0];
	p_j = coords[1];
	p_k = coords[2];
	
	libsci_acc_HostAlloc((void**)&A, sizeof(val_type)*N*N);
	libsci_acc_HostAlloc((void**)&B, sizeof(val_type)*N*N);
	libsci_acc_HostAlloc((void**)&B, sizeof(val_type)*N*N);
	
	
	#ifdef DEBUGOUTPUT 
		std::cout << "Hello from Rank " << rank << " with coordinates (i,j,k) = (" << p_i <<","<<p_j<<","<<p_k<<")\n";
	#endif
	
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
