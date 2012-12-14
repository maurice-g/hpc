#include <iostream>
#include "DistGEMM.hpp"
#include <cstdlib>
#include <omp.h>
int main(int argc, char* argv[]) {
	int N = std::atoi(argv[1]);
	int nprocs = std::atoi(argv[2]);
	int topsize = std::atoi(argv[3]);	//size of 3dtopology: topsize x topsize x topsize
	DistGEMM m(N,nprocs,topsize);

	return 0;
}
