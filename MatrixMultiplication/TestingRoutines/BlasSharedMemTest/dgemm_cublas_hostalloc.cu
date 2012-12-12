//dgemm testing cublas
//Stefan H
//5.12.12
#include <cublas.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <matrix.hpp>
#include <fillMatrix.hpp>
#include <Timer.hpp>
#include <Measurement.hpp>

	

int main() {

	typedef hpc12::matrix<double,hpc12::column_major> matrix_type;
	for (int N = 512;N < 10000;N*=2) {
		matrix_type A(N,N);
		matrix_type B(N,N);
		matrix_type C(N,N);
		
		double * d_A, * d_B, *d_C;
		cudaMalloc((void**) &d_A, N*N*sizeof(double));
                cudaMalloc((void**) &d_B, N*N*sizeof(double));
                cudaMalloc((void**) &d_C, N*N*sizeof(double));

		fillMatrix(A);
		fillMatrix(B);

		//see whether it is faster to use pinned (page locked) memory for matrices on host
		double * h_A, *h_B, *h_C;
		cudaMallocHost((void**) &h_A,N*N*sizeof(double));
		cudaMallocHost((void**) &h_B,N*N*sizeof(double));
		cudaMallocHost((void**) &h_C,N*N*sizeof(double));


		//transfer data into pinned memory
		cudaMemcpy(h_A, A.data(), N*N*sizeof(double),cudaMemcpyHostToHost);
		cudaMemcpy(h_B, B.data(), N*N*sizeof(double),cudaMemcpyHostToHost);

		// include time to copy to /from device
		Timer _t(1);

		cudaMemcpy (d_A, h_A,N*N*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy (d_B, h_B ,N*N*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemset (d_C, 0., N*N*sizeof(double));
		char trans = 'N';
		double alpha = 1.;
		double beta = 0.;
		cublasDgemm(trans,trans,N,N,N,alpha,d_A,N,d_B,N,beta,d_C,N);
		cudaMemcpy(h_C,d_C,N*N*sizeof(double),cudaMemcpyDeviceToHost);
		
		_t.stop();
		Measurement m("cublasDgemm_MallocHost",N,N,_t.elapsed_s());
		std::cout << m;
		//free memory on device
		cudaFree(d_A);
		cudaFree(d_B);
		cudaFree(d_B);

		//free pinned memory
		cudaFreeHost(h_A);cudaFreeHost(h_B);cudaFreeHost(h_C);
	}


	return 0;
}
