//dgemm testing cublas
//Stefan H
//5.12.12
#include <cublas.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <matrix.hpp>
#include <fillMatrix.hpp>
#include <Timer.hpp>
#include <chrono>

	
void
dgemm_cublas(MatrixType &A, MatrixType &B, MatrixType & C) {
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
        cublasDgemm(transA,transA,M,N,K,alpha,A.data(),LDA,B.data(),LDB,beta,C.data(),LDC);

}


int main() {

	typedef hpc12::matrix<double,hpc12::column_major> matrix_type;
	for (int N = 512;N < 20000;N*=2) {
		matrix_type A(N,N);
		matrix_type B(N,N);
		matrix_type C(N,N);
		fillMatrix(A);
		fillMatrix(B);
		//C = A*B
		Timer _t(1);
		dgemm_libsci(A,B,C);
		_t.stop();
		Measurement m("dgemm with libsci,?noacc?,16thrds",N,N,_t.elapsed_s());
		std::cout << m;
	}


	return 0;
}
