/*!
* \author Stefan H.
* \date 26.11.2012
* HPC CSE project
*/

#ifndef MEAS_HPP
#define MEAS_HPP
#include <string>
#include <iostream>

/*
* \class Measurement
* \brief Class implementing the results of a measurement
*/
class Measurement {
	public:
		/*!
		* Constructor
		* @param description is a string which holds a description of the measurement
		* @param m is the number of rows of the matrix
		* @param n is the number of columns of the matrix
		* @param runtime
		*/
		Measurement(std::string description, std::size_t m, std::size_t n, double runtime) :
					description_(description),
					m_(m),
					n_(n),
					runtime_(runtime) {
					flops_ = 2*m*n*m / runtime_;	//change this to the exact value
					};
	
	
	private:
		std::string description_;
		std::size_t m_,n_;
		double runtime_;
		double flops_;
		
	friend std::ostream & operator << (std::ostream& os, const Measurement& x);



};

/*!
* Function which allows us to write std::ostream << measurement and write the result directly into a file
* @param os : ostream object
* @param x : a Measurement object
*/
std::ostream& 
operator << (std::ostream& os, const Measurement& x) {
	os << x.description_ << " " << x.m_ << " " << x.n_ << " " << x.runtime_  << " " << x.flops_ << "\n";
	return os;
}
#endif
