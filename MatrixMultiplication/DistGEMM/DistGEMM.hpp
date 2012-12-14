#ifndef MATRIXMULT
#define MATRIXMULT

/*
* \brief Class implementing distributed gemm for a n x n matrix
*/

class DistGEMM {
public:
	typedef double 		val_type;
	typedef unsigned int	count_type;
	#DEFINE mpi_val_type	MPI_DOUBLE

	/*!
	* Constructor
	* @param n is the size of the entire matrix
	*/
	DistGEMM(int N, int P);					// constructor: initialize communicators, matrix allocation
	initializeLehmer();					// initialize Lehmer matrix A (where k=0) and B (where i=0)
	performGEMM();						// perform the matrix multiplication

private:
	MPI_Comm	comm_i, comm_j, comm_k;
	val_type *A, *B, *C;

	count_type blocksize;					// matrix size per node
	count_type P;						// # of nodes
	count_type p_i,p_j,p_k;					// topology index of p	

};


#endif
