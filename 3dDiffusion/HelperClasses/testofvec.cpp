
#include <iostream>
#include "Vector3d.hpp"
int main() {
	unsigned int N = 4;
	Vector3d<double> V(N,N,N);
	
	V(1,1,1) = 1;
	std::cout << V.get_sizeX()<< "\n";
	std::cout << V;

}