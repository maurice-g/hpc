//dgemm testing cublas, measuring flops with and without copy to dev and membdwdth
//Stefan H
//11.12.12
#include <cublas.h>
#include <cuda_runtime_api.h>
#include <iostream>
#include <matrix.hpp>
#include <fillMatrix.hpp>
#include <Timer.hpp>
#include <Measurement.hpp>

#ifndef FLOPSWITHOUTHCPY
#define FLOPSWITHOUTCPY 0
#endif

#ifndef MEMBDWDTH
#define MEMBDWDTH 0
#endif

#ifndef FLOPS
#define FLOPS 0
#endif
	

int main() {
	std::cout << "performance  of cublas dgemm using pinned memory \n";
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
		
		cudaEvent_t start,stop;
		cudaEventCreate(&start);
		cudaEventCreate(&stop);
		//cudaEventRecord(start,0);

		cudaMemcpy (d_A, h_A,N*N*sizeof(double),cudaMemcpyHostToDevice);
		cudaMemcpy (d_B, h_B ,N*N*sizeof(double),cudaMemcpyHostToDevice);

		//cudaEventRecord(stop,0);
		//cudaEventSynchronize(stop);
		float elapsed_time;
		double elapsed;
		/*cudaEventElapsedTime(&elapsed_time,start,stop);//time in miliseconds, bw in bytes: -->GB/s: 1e-9*1e3
		std::cout << N  << "Peak memory bandwidth H2D(GB/s): " << 1e-6*2*N*N*sizeof(double)/elapsed_time << "\n";
		cudaEventDestroy(start);
		cudaEventDestroy(stop);
*/
		cudaMemset (d_C, 0., N*N*sizeof(double));
		char trans = 'N';
		double alpha = 1.;
		double beta = 0.;
		cudaEventRecord(start,0);
		cublasDgemm(trans,trans,N,N,N,alpha,d_A,N,d_B,N,beta,d_C,N);
		cudaEventRecord(stop,0);
		cudaEventSynchronize(stop);
		cudaEventElapsedTime(&elapsed_time,start,stop);elapsed=elapsed_time;
		
		//craaazy bug: try double(N*N*N)/(elapsed*0.001) ---> result is 0 WTF??
		std::cout << N << " runtime: "<< elapsed  << " peak FLOPS without copying data: " << double(N)*double(N)*double(N)/(elapsed*0.001) <<"\n";
		

		cudaMemcpy(h_C,d_C,N*N*sizeof(double),cudaMemcpyDeviceToHost);		

		//free memory on device
		cudaFree(d_A);
		cudaFree(d_B);
		cudaFree(d_B);

		//free pinned memory
		cudaFreeHost(h_A);cudaFreeHost(h_B);cudaFreeHost(h_C);
	}


	return 0;
}
