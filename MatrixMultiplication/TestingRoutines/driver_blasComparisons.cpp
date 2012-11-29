#include <iostream>
#include "HelperClasses/matrix.hpp"
#include "HelperClasses/dgemm_libsci.hpp"
#include "HelperClasses/fillMatrix.hpp"
#include "HelperClasses/Measurement.hpp"

#include <chrono>
#include <omp.h>

int main() {
	omp_set_num_threads(16);
	typedef hpc12::matrix<double,hpc12::column_major> matrix_type;
	for (int N = 512;N < 10000;N*=2) {
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
