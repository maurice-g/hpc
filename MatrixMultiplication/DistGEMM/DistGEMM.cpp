#include "DistGEMM.hpp"
#include <mpi.h>
#include <libsci_acc.h>
#include <cassert>
#include <cmath>
#include <iostream>

DistGEMM::DistGEMM(int n, int numprocs, int cubes) {
	assert(cubes*cubes*cubes == numprocs);	
	assert(n % cubes == 0);
	
	// removed size,rank from class, isn't size==numprocs??
	int size;	
	MPI_Comm_size(MPI_COMM_WORLD,&size);
	assert(size==numprocs);

	libsci_acc_init();

	P = numprocs;
	cubeSize = cubes;					//# processors in each direction of cartesian topology
	blocksize = n/cubeSize;
	
	//build MPI topology
	int nums[3] = {cubes, cubes, cubes};			//number of nodes in each dimension (x,y,z)
	MPI_Dims_create(size,3,nums);
	//build global communicator, then subdivide into comm_i, comm_j, comm_k
	int periodic[3] = {0, 0, 0};
	MPI_Comm cart_comm;
	MPI_Cart_create(MPI_COMM_WORLD,3,nums,periodic,true,&cart_comm);
	
	//assign comm_i, comm_j, comm_k
	int logv_i[] = {true,false,false};
	int logv_j[] = {false,true,false};
	int logv_k[] = {false,false,true};
	MPI_Cart_sub(cart_comm,logv_i,&comm_i);
	MPI_Cart_sub(cart_comm,logv_j,&comm_j);
	MPI_Cart_sub(cart_comm,logv_k,&comm_k);

        MPI_Comm_rank(MPI_COMM_WORLD,&rank_cart);

	//get and save the new ranks for each communications layer
	MPI_Comm_rank(comm_i,&rank_i);
	MPI_Comm_rank(comm_j,&rank_j);
	MPI_Comm_rank(comm_k,&rank_k);
	//std::cout << "Global= " << rank << ";i=" << ranki << ";j=" << rankj << ";k=" << rankk << std::endl;

	//assign the coordinates to each rank
	int d = 3;
	int coords[3];
	MPI_Cart_coords(cart_comm,rank_cart,d,coords);
	p_i = coords[0];
	p_j = coords[1];
	p_k = coords[2];
	
	// get the root of each 1D layer
	int rootcoord = 0;
	// ATTENTION: root_i!=p_i(0) 
	MPI_Cart_rank(comm_i, &rootcoord, &root_i);
	MPI_Cart_rank(comm_j, &rootcoord, &root_j);
	MPI_Cart_rank(comm_k, &rootcoord, &root_k);

	libsci_acc_HostAlloc((void**)&A, sizeof(val_type)*blocksize*blocksize);
	libsci_acc_HostAlloc((void**)&B, sizeof(val_type)*blocksize*blocksize);
	libsci_acc_HostAlloc((void**)&C, sizeof(val_type)*blocksize*blocksize);
	
	#ifdef DEBUGOUTPUT 
		std::cout << "Hello from Rank " << rank_cart << " with coordinates (i,j,k) = (" << p_i <<","<<p_j<<","<<p_k<<")\n";
	#endif
	
}

void DistGEMM::initializeLehmer() {
	// fill A and B as Lehmer matrices
	for (count_type j=0; j<blocksize; j++) {
		for (count_type i=0; i<blocksize; i++) {
			A[j*blocksize+i] = val_type(std::min(p_i*blocksize+i+1,p_j*blocksize+j+1))/std::max(p_i*blocksize+i+1,p_j*blocksize+j+1);
			B[j*blocksize+i] = val_type(std::min(p_j*blocksize+i+1,p_k*blocksize+j+1))/std::max(p_j*blocksize+i+1,p_k*blocksize+j+1);
		}
	}
}

void DistGEMM::performGEMM() {
	// send root_k's matrix A to all other processors in communicator comm_k
	MPI_Bcast(A, blocksize*blocksize, mpi_val_type, root_k, comm_k);
	// send root_i's Matrix B to all other processors in communicator comm_i
	MPI_Bcast(B, blocksize*blocksize, mpi_val_type, root_i, comm_i);
	
	char transA = 'N';
	double alpha = 1.;
	double beta = 0.;
	
	dgemm(transA, transA, blocksize, blocksize, blocksize, alpha, A, blocksize, B, blocksize, beta, C, blocksize);
	MPI_Reduce((rank_j==root_j ? MPI_IN_PLACE : C), C, blocksize*blocksize, mpi_val_type, MPI_SUM, root_j, comm_j);
}

void DistGEMM::output_result() {
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0; i<cubeSize; i++) {
		for (int k=0; k<cubeSize; k++) {
			MPI_Barrier(MPI_COMM_WORLD);
			// very strange output behaviour: always (0,0) first, other procs in random order
			if (p_i==i && p_k==k && rank_j==root_j) {
				std::cout << "output_from(i,k): " << p_i << "," << p_k << std::endl;
				for (int op=0; op<blocksize*blocksize; op++) {
					std::cout << C[op] << std::endl;
				}
			}
		}
	}
}
			

DistGEMM::~DistGEMM() {
	libsci_acc_FreeHost(A);
	libsci_acc_FreeHost(B);
	libsci_acc_FreeHost(C);
}	
