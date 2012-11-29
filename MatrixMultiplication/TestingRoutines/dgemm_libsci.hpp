// dgemm version calling libsci library (NOT accelerated)
//Stefan
//C = A*B


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
