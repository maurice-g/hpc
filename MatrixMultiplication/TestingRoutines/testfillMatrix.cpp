//Just a small test for the fillMatrix function
//Stefan, 28.11.12

#include <iostream>
#include "../HelperClasses/fillMatrix.hpp"
#include "../HelperClasses/matrix.hpp"



int main() {
	typedef hpc12::matrix<double,hpc12::column_major> matrix_type;
	matrix_type A(10,10);
	fillMatrix(A);
	
	std::cout << A;

	return 0;
}
