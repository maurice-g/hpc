
#ifndef MATRIXMULT
#define MATRIXMULT

#include <mpi.h>
#include <string>

/*
* \class DistGEMM
* \brief Class implementing distributed gemm for a n x n matrix
*/
/**
* This class implements a distributed (MPI) gemm using the 3d-matrix-multiplicaiton
* algorithm described in....
* there are several possibilities for setting up the matrices:
* \li load matrices from file, \sa setup
* \li -generate matrices on each node individually (lehmer matrices)
*/


class DistGEMM {
public:
	typedef double 		val_type;
	typedef unsigned int	count_type;
	#define mpi_val_type	MPI_DOUBLE

	/*!
	* Constructor, ! Requirements: N mod d == 0, d^3 == P
	* @param N is the size of the matrix
	* @param P is the number of nodes
	* @param d is the size of the 3-dimensional node-topology
	* 
	*/
	DistGEMM(int N, int P, int d);					// constructor: initialize communicators & topology, matrix allocation
	
	/*! 
	* Destructor: Frees Memory of matrices A,B and C
	*/
	~DistGEMM();
	
	
	/*!
	* testing function \n
	* each node fills its (local) matrix part in a way
	* such that the global matrices are lehmer matrices
	*/ 
	void initializeLehmer();					// initialize Lehmer matrix A (where k=0) and B (where i=0)
	
	
	
	/*!
	* This function performs C = A*B
	* using the distributed 3d-matrix-multiplication algorithm
	* the local computation on each node are performed using the 
	* libsci_acc version of dgemm (which distributes work between gpu-cpu)
	*/
	void performGEMM();						// perform the matrix multiplication
	
	/*!
	*This function writes the matrix C (the result of A*B) to a file using the following format:\n
	* Column major storage, but as a single array (column)
	*@param filename
	**/
	void output_result(std::string filename);
	
		
	/*!
	* Use this function to load the matrices A,B from files (filenameA,filenameB)\n
	* Required Format: plain-text files, NxN entries in one column\n
	* these numbers are interpreted as the matrix entries stored in column major format
	**/
	void setup(std::string filenameA, std::string filenameB);						

private:
	MPI_Comm cart_comm, comm_i, comm_j, comm_k;
	int rank_cart,rank_i,rank_j,rank_k;				// the rank of each communicator
	val_type *A, *B, *C;

	count_type blocksize;					// matrix size per node
	count_type P;						// # of nodes
	count_type p_i,p_j,p_k;					// topology index of p	
	count_type cubeSize;					//MPI 3d Topology length ->(P^1/3)
	int root_i, root_j, root_k;

	//function gets called from setup()
	void readMatrix(std::string filename, val_type *matrix);
};


#endif
