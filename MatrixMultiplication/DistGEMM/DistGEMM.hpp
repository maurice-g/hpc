#ifndef MATRIXMULT
#define MATRIXMULT

#include <mpi.h>

/*
* \brief Class implementing distributed gemm for a n x n matrix
*/

class DistGEMM {
public:
	typedef double 		val_type;
	typedef unsigned int	count_type;
	#define mpi_val_type	MPI_DOUBLE

	/*!
	* Constructor
	* @param N is the size of the entire matrix
	* @param P is the number of nodes 
	*/

	DistGEMM(int N, int P);					// constructor: initialize communicators & topology, matrix allocation
	void initializeLehmer();					// initialize Lehmer matrix A (where k=0) and B (where i=0)
	void performGEMM();						// perform the matrix multiplication


private:
	int rank, size;
	MPI_Comm	comm_i, comm_j, comm_k;
	val_type *A, *B, *C;

	count_type blocksize;					// matrix size per node
	count_type P;							// # of nodes
	count_type p_i,p_j,p_k;					// topology index of p	
	count_type cubeSize;					//MPI 3d Topology length ->(P^1/3)
	count_type root_i, root_j, root_k;

};


#endif
