#include "Diffusion3D.hpp"
#include <mpi.h>
#include <cassert>
#include <iostream>

//constructor
Diffusion3D::Diffusion3D(val_type dx, count_type nx, count_type ny, count_type nz,coord_type topology,val_type D, val_type T):
						D_(D),T_(T),dx_(dx),
						global_nx_(nx),global_ny_(ny),global_nz_(nz),
						global_lengthx_(nx*dx),global_lengthy_(ny*dx),global_lengthz_(nz*dx)
{
	
	std::cout << "Constructor Called\n";
	//check # nodes vs mesh points: each node has the same number of meshpoints!
	assert(nx % topology[0] == 0);
	assert(ny % topology[1] == 0);
	assert(nz % topology[2] == 0);
	
	
	
	dt_ = dx_*dx_/(6*D);	// dt < dx^2 /(4*D) for stability of FTCS

	density_.resize(2,2,2,0);
}


//Destructor
Diffusion3D::~Diffusion3D() {
	//free mpi types!
	std::cout << "Destructor Called\n";
	
}


//function which writes debug information to os
void Diffusion3D::write_debug_info(std::ostream &os) const {
	os << "\n-------------Info on object at (this-pointer) " << this <<"------------\n";
	os << "SimParameters: \t dx:\t\t\t" << this->dx_ << "\n"
		<< "\t\t dt:\t\t\t" << this->dt_ << "\n"
		<< "\t\t D:\t\t\t"  << this->D_ << "\n"
		<< "\t\t T:\t\t\t"  << this->T_ << "\n"
		<< "\t\t globallengthx:\t\t" << this->global_lengthx_ << "\n"
		<< "\t\t globallengthy:\t\t" << this->global_lengthy_ << "\n"
		<< "\t\t globallengthz:\t\t" << this->global_lengthz_ << "\n"
		<< "\t\t global nx: \t\t"  << this->global_nx_ << "\n"
		<< "\t\t global ny: \t\t"  << this->global_ny_ << "\n"
		<< "\t\t global nz: \t\t"  << this->global_nz_ << "\n\n";
		
	os << "MPIParameters: \t Size: \t\t\t" << this->size_ <<"\n"
		<< "\t\t Rank \t\t\t" << this->rank_ <<"\n"
		<< "\t\t CartCoords(x,y,z):\t"<<this->cartesian_coords_[0] << " " << this-> cartesian_coords_[1] << " " << this->cartesian_coords_[2]<<"\n"
		<< "\t\t local nx:\t\t" << this->local_nx_ << "\n"
		<< "\t\t local ny:\t\t" << this->local_ny_ << "\n"
		<< "\t\t local nz:\t\t" << this->local_nz_ << "\n"
		<< "\t\t left neighbour:\t"<<this->left_ <<"\n"
		<< "\t\t right neighbour: \t"<<this->right_<<"\n"
		<< "\t\t front neighbour:\t"<<this->front_ <<"\n"
		<< "\t\t back neighbour: \t"<<this->back_<<"\n"
		<< "\t\t top neighbour:\t\t"<<this->top_ <<"\n"
		<< "\t\t bottom neighbour: \t"<<this->bottom_<<"\n\n"
		;
		
		
	os << "Orientation of coord system to understand left,right,top etc, followed by densities\n" << this->density_ <<"\n";
	os << "Check whether density_ vector is aligned to 64 bytes:\n"
		<< "\t\t Starting Adress\t\t" << &this->density_(0,0,0) << "\n"
		<< "\t\t Adress of (1,0,0)\t\t" << &this->density_(1,0,0) <<"\n"
		<< "\t\t Adress of (0,1,0)\t\t" << &this->density_(0,1,0) <<"\t Stride in y-direction:\t"<<int(&this->density_(0,1,0))-int(&this->density_(0,0,0))<<" Bytes\n"
		<< "\t\t Adress of (0,0,1)\t\t" << &this->density_(0,0,1) <<"\t Stride in z-direction:\t"<<int(&this->density_(0,0,1))-int(&this->density_(0,0,0))<<" Bytes\n";
	os << "Stride in y-direction should be:\t\t" << density_.get_sizeX()*sizeof(val_type) <<"\n";
	os << "Stried in z-direction should be:\t\t" << density_.get_sizeX()*density_.get_sizeY()*sizeof(val_type)<<"\n";
	os << "\n------------------------------------------------------------------------\n\n";
}


