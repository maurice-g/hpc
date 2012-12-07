#include <iostream>
#include <matrix.hpp>
#include <fillMatrix.hpp>
#include <Measurement.hpp>

#include <chrono>
#include <omp.h>



extern "C" void dgemm_ (char & transa, char & transb,
                        int & m, int & n, int & k,
                        double & alpha, double * A, int & LDA,
                        double * B, int & LDB,
                        double & beta, double * C, int & LDC);

template<class MatrixType>
void
dgemm_libsci(MatrixType &A, MatrixType &B, MatrixType & C) {
        //assert A columnmajor
        char transA = 'N';
        int M = A.num_rows();
        int N = B.num_cols();
        int K = A.num_cols();
        double alpha = 1.;
        int LDA = A.leading_dimension();
        int LDB = B.leading_dimension();
        double beta = 0.;
        int LDC = C.leading_dimension();
        dgemm_(transA,transA,M,N,K,alpha,A.data(),LDA,B.data(),LDB,beta,C.data(),LDC);


}



int main() {
	omp_set_num_threads(16);
	typedef hpc12::matrix<double,hpc12::column_major> matrix_type;
	for (int N = 512;N < 20000;N*=2) {
		matrix_type A(N,N);
		matrix_type B(N,N);
		matrix_type C(N,N);
		fillMatrix(A);
		fillMatrix(B);
		//C = A*B
		std::chrono::time_point<std::chrono::high_resolution_clock> start,end;
		start = std::chrono::high_resolution_clock::now();
		dgemm_libsci(A,B,C);
		end = std::chrono::high_resolution_clock::now();
		double elapsed_seconds = std::chrono::duration<double>(end-start).count();
		Measurement m("dgemm with libsci,?noacc?,16thrds",N,N,elapsed_seconds);
		std::cout << m;
	}


	return 0;
}
