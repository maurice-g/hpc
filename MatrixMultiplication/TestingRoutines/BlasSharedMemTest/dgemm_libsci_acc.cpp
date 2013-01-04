#include <iostream>
#include <matrix.hpp>
#include <fillMatrix.hpp>
#include <Measurement.hpp>

#include <Timer.hpp>
#include <libsci_acc.h>
#include <omp.h>

int main() {
	omp_set_num_threads(16);
	libsci_acc_init();
	double *A, *B, *C;
	char transA = 'N';
	double alpha = 1.;
	double beta = 0.;
	for (int N = 512;N < 25000;N+=512) {
                libsci_acc_HostAlloc( (void **)&A, sizeof(double)* N*N  );
                libsci_acc_HostAlloc( (void **)&B, sizeof(double)* N*N );
                libsci_acc_HostAlloc( (void **)&C, sizeof(double)* N*N  );
		//fillMatrix(A);
		for (int i=0; i<N; i++) {
			for (int j=0; j<N; j++) {
				A[i*N+j] = i*j;
				B[i*N+j] = i*j;
			}
		}
		//fillMatrix(B);
		//C = A*B
		Timer t(1);
		dgemm( transA, transA, N, N, N, alpha, A, N, B, N, beta, C,N );
		t.stop();
		Measurement m("dgemm with libsci_acc_16thrds",N,N,t.elapsed_s());
		std::cout << m;
      		libsci_acc_FreeHost( A );
      		libsci_acc_FreeHost( B );
      		libsci_acc_FreeHost( C );
	}


	return 0;
}
