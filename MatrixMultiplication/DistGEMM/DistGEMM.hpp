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
	* @param N is the size of the entire matrix
	* @param P is the number of nodes 
	*/
<<<<<<< HEAD
	DistGEMM(int N, int P);					// constructor: initialize communicators, matrix allocation
	void initializeLehmer();					// initialize Lehmer matrix A (where k=0) and B (where i=0)
	void performGEMM();						// perform the matrix multiplication
=======
	DistGEMM(int N, int P);					// constructor: initialize communicators & topology, matrix allocation
	initializeLehmer();					// initialize Lehmer matrix A (where k=0) and B (where i=0)
	performGEMM();						// perform the matrix multiplication
>>>>>>> 6883e97ec3b4dc58db2292bc3782afbe4495f1c2

private:
	MPI_Comm	comm_i, comm_j, comm_k;
	val_type *A, *B, *C;

	count_type blocksize;					// matrix size per node
	count_type P;						// # of nodes
	count_type p_i,p_j,p_k;					// topology index of p	

};


#endif
