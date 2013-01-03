#include "DistGEMM.hpp"
#include <mpi.h>
#include <vector>
#include <libsci_acc.h>
#include <cassert>
#include <cmath>
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
	MPI_Cart_create(MPI_COMM_WORLD,3,nums,periodic,true,&cart_comm);
	
	//assign comm_i, comm_j, comm_k
	int logv_i[] = {true,false,false};
	int logv_j[] = {false,true,false};
	int logv_k[] = {false,false,true};
	MPI_Cart_sub(cart_comm,logv_i,&comm_i);
	MPI_Cart_sub(cart_comm,logv_j,&comm_j);
	MPI_Cart_sub(cart_comm,logv_k,&comm_k);

        MPI_Comm_rank(cart_comm,&rank_cart);

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

void DistGEMM::output_result(std::string filename) {
	MPI_Comm comm_gather;
	int dims[] = {true,false,true};
	int my_rank;
	
	MPI_Cart_sub(cart_comm,dims,&comm_gather);
	MPI_Comm_rank(comm_gather, &my_rank);
	
	if (p_j!=0)
		return;
	
	val_type *result=0x0;
	if (my_rank==0) // the proc receiving all data
		result=new val_type[blocksize*blocksize*cubeSize*cubeSize];
	
	std::vector<int> recvcounts(cubeSize*cubeSize, blocksize*blocksize);
	std::vector<int> displs(cubeSize*cubeSize, 0);
	for (count_type i=0; i<cubeSize; i++) {
		for (count_type k=0; k<cubeSize; k++) {
			int coords[] = {i,k};
			int rank_gather; 
			MPI_Cart_rank(comm_gather, coords, &rank_gather);
			displs[rank_gather] = (blocksize*blocksize)*(i+cubeSize*k); //store blocks column major
		}
	}
	
	MPI_Gatherv(C, blocksize*blocksize, mpi_val_type, result, &recvcounts[0], &displs[0], mpi_val_type, 0, comm_gather);
	
	if (my_rank==0) {
		std::ofstream outfile;
		std::stringstream matSizeInfo;
		matSizeInfo << blocksize*cubeSize;
		filename += matSizeInfo.str();
		filename += ".txt";
		outfile.open(filename.c_str(),std::ios_base::app | std::ios_base::out);
		if (outfile.is_open()) {
			for (int k=0; k<cubeSize; k++) {
				for (int column=0; column<blocksize; column++) {
					for (int i=0; i<cubeSize; i++) {
						for (int op=0; op<blocksize; op++) {
							outfile << result[(blocksize*blocksize)*(cubeSize*k+i)+column*blocksize+op] << "\n";
						}
					}
				}
			}
		} else {
			std::cout << "Error: file not open\n";
		}
		outfile.close();
	}
	delete[] result;
}

void DistGEMM::readMatrix(std::string filename, val_type *matrix) {
	// reads and stores matrix blockwise into large array (column major)
	count_type N=blocksize*blocksize*cubeSize*cubeSize;
	val_type line;
	
	count_type iblock=0;
	count_type jblock=0;
	count_type column=0;
	count_type op=0;
	count_type c=0;	

	std::ifstream inputfile(filename.c_str());
	if (inputfile.is_open()) {
		while(inputfile.good()) {
			c++;
			
			inputfile >> line;
			matrix[blocksize*blocksize*(iblock+cubeSize*jblock)+blocksize*column+op] = line;
                        
			if (++op==blocksize) {
                                op=0;
                                if (++iblock==cubeSize) {
                                        iblock=0;
                                        if (++column==blocksize) {
                                                column=0;
                                                if (++jblock==cubeSize) break;
                                        }
                                }
                        }
		}
	}

	else std::cerr << "Error while reading file " << filename << std::endl;
	if (c != N) std::cerr << "Warning: number of matrix elements given to class (" << N << ") is not equal to number of matrix elements in file " << filename << "(" << c << ")" << std::endl;
	inputfile.close();
}

void DistGEMM::setup(std::string filenameA, std::string filenameB) {
	count_type N = blocksize*blocksize*cubeSize*cubeSize;
	MPI_Request *reqs_send_A=new MPI_Request[cubeSize*cubeSize];
	MPI_Status *status_send_A=new MPI_Status[cubeSize*cubeSize];
        MPI_Request *reqs_send_B=new MPI_Request[cubeSize*cubeSize];
	MPI_Status *status_send_B=new MPI_Status[cubeSize*cubeSize];
	val_type *matrixA, *matrixB;
	matrixA=0x0;
	matrixB=0x0;
	
	if (rank_cart==0) {
		// reads matrix A 
		matrixA = new val_type[N];
		readMatrix(filenameA, matrixA);
		
		for (count_type i=0; i<cubeSize; i++) {
			for (count_type j=0; j<cubeSize; j++) {
				int coords[] = {i,j,0};
				int sendTo;
				MPI_Cart_rank(cart_comm, coords, &sendTo);
				MPI_Isend(&matrixA[(blocksize*blocksize)*(cubeSize*j+i)], blocksize*blocksize, mpi_val_type, sendTo, 0, cart_comm, &reqs_send_A[i*cubeSize+j]);
			}
		}
	}
	
	if (rank_cart==1) {
                // reads matrix B
                matrixB = new val_type[N];
                readMatrix(filenameB, matrixB);

                for (count_type j=0; j<cubeSize; j++) {
                        for (count_type k=0; k<cubeSize; k++) {
                                int coords[] = {0,j,k};
                                int sendTo;
                                MPI_Cart_rank(cart_comm, coords, &sendTo);
                                MPI_Isend(&matrixB[(blocksize*blocksize)*(cubeSize*k+j)], blocksize*blocksize, mpi_val_type, sendTo, 0, cart_comm, &reqs_send_B[j*cubeSize+k]);
                        }
                }
        }
	
	MPI_Request reqs_recv_A, reqs_recv_B;
	MPI_Status status_recv_A, status_recv_B;

	if (p_k==0) MPI_Irecv(A, blocksize*blocksize, mpi_val_type, 0, 0, cart_comm, &reqs_recv_A);
	if (p_i==0) MPI_Irecv(B, blocksize*blocksize, mpi_val_type, 1, 0, cart_comm, &reqs_recv_B);

	if (rank_cart==0) MPI_Waitall(cubeSize*cubeSize, &reqs_send_A[0], &status_send_A[0]);
	if (rank_cart==1) MPI_Waitall(cubeSize*cubeSize, &reqs_send_B[0], &status_send_B[0]);
	if (p_k==0) MPI_Waitall(1, &reqs_recv_A, &status_recv_A);
	if (p_i==0) MPI_Waitall(1, &reqs_recv_B, &status_recv_B);

	delete[] reqs_send_A, status_send_A, reqs_send_B, status_send_B, matrixA, matrixB;
	
	MPI_Barrier(MPI_COMM_WORLD);
}

/*
void DistGEMM::output_result() {
	MPI_Barrier(MPI_COMM_WORLD);
	for (int i=0; i<cubeSize; i++) {
		MPI_Barrier(MPI_COMM_WORLD);
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
*/			

DistGEMM::~DistGEMM() {
	libsci_acc_FreeHost(A);
	libsci_acc_FreeHost(B);
	libsci_acc_FreeHost(C);
}	
