#include "Diffusion3D.hpp"
#include <mpi.h>
#include <cassert>
#include <iostream>

//constructor
Diffusion3D::Diffusion3D(val_type dx, count_type nx, count_type ny, count_type nz, coord_type a, coord_type b,coord_type topology,val_type D, val_type T):
						D_(D),T_(T),dx_(dx),
						global_nx_(nx),global_ny_(ny),global_nz_(nz)
{
	
	std::cout << "Constructor Called\n";
	//check # nodes vs mesh points: each node has the same number of meshpoints!
	assert(nx % topology[0] == 0);
	assert(ny % topology[1] == 0);
	assert(nz % topology[2] == 0);
	
	//check for domain: a is front left down, b is back right up if x right, y up, z front
	assert(a[0] < b[0]);
	assert(a[1] < b[1]);
	assert(a[2] > b[2]);
	
	dt_ = dx_*dx_/(6*D);
	


}

Diffusion3D::~Diffusion3D() {
	//free mpi types!
	std::cout << "Destructor Called\n";
	
}
