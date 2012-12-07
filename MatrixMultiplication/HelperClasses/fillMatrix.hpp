/*!
HPCSE project dgemm
* \author Stefan H
* \date 28.11.2012
* 
*/

#ifndef FILLMAT_HPP
#define FILLMAT_HPP

#include <cstring>

/*! 
* \p Function which fills a matrix: A(i,j) = i*j
* - Concepts for Type MatrixType:
* - A.num_cols() exists and returns the number of columns of A convertible to size_t
* - A.num_rows() exists and returns the number of rows of A convertible to size_t
* - operator(i,j) is overloaded and gives access to the underlying element A(i,j) (rhs assignable)
* @param &A: Matrix which you want to fill
*
*/

template<class MatrixType>
void
fillMatrix(MatrixType &A) {
	
		for (std::size_t j = 0; j < A.num_cols(); j++) {
			for (std::size_t i = 0; i < A.num_rows(); i++) {
				A(i,j) = i*j;
			}
		}
}



#endif
