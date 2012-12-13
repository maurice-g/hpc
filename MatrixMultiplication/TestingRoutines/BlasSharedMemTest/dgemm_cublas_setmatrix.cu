//dgemm testing cublas version 2, using set /get matrix for transferring data d2h h2d
//Stefan H
//13.12.12
#include <cublas_v2.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <matrix.hpp>
#include <fillMatrix.hpp>
#include <Timer.hpp>
#include <Measurement.hpp>

	

int main() {
	cublasHandle_t handle;
	cublasStatus_t status = cublasCreate(&handle);

	typedef hpc12::matrix<double,hpc12::column_major> matrix_type;
	for (int N = 512;N < 10000;N*=2) {
		matrix_type A(N,N);
		matrix_type B(N,N);
		matrix_type C(N,N);
		
		fillMatrix(A);
		fillMatrix(B);

		
		double * d_A, * d_B, *d_C;
		cudaMalloc((void**) &d_A, N*N*sizeof(double));
        cudaMalloc((void**) &d_B, N*N*sizeof(double));
        cudaMalloc((void**) &d_C, N*N*sizeof(double));

		//allocate pinned memory
		double * h_A, *h_B, *h_C;
		cudaMallocHost((void**) &h_A,N*N*sizeof(double));
		cudaMallocHost((void**) &h_B,N*N*sizeof(double));
		cudaMallocHost((void**) &h_C,N*N*sizeof(double));



		//transfer data into pinned memory
		cudaMemcpy(h_A, A.data(), N*N*sizeof(double),cudaMemcpyHostToHost);
		cudaMemcpy(h_B, B.data(), N*N*sizeof(double),cudaMemcpyHostToHost);

		
		Timer _t(1);

		cublasSetMatrix (A.num_rows(),A.num_cols(),sizeof(double),h_A,A.leading_dimension(),d_A,A.leading_dimension());
		cublasSetMatrix (B.num_rows(),B.num_cols(),sizeof(double),h_B,B.leading_dimension(),d_B,B.leading_dimension());
		cudaMemset (d_C, 0., N*N*sizeof(double));
		//C = A*B
		double alpha = 1.;
		double beta = 0.;
		
		cublasDgemm(handle,CUBLAS_OP_N,CUBLAS_OP_N,N,N,N,&alpha,d_A,N,d_B,N,&beta,d_C,N);
		
		cublasGetMatrix(C.num_rows(),C.num_cols(),sizeof(double),d_C,C.leading_dimension(),h_C,C.leading_dimension());
		
		
		_t.stop();
		Measurement m("cublasDgemm v2, setmatrix etc",N,N,_t.elapsed_s());
		std::cout << m;
		//free memory on device
		cudaFree(d_A);
		cudaFree(d_B);
		cudaFree(d_C);
		
		cudaFreeHost(h_A);cudaFreeHost(h_B);cudaFreeHost(h_C);
	}
	cublasDestroy(handle);

	return 0;
}
