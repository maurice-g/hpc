
#include <iostream>
#include "Vector3d.hpp"
int main() {
	unsigned int N = 4;
	Vector3d<double> V(N,N,N);
	
	V(1,1,1) = 1;
	V(1,2,3) = 123;
	V(0,0,0) = -1;
	std::cout << V.get_sizeX()<< "\n";
	std::cout << V;

}
