#include "Diffusion3D.hpp"
#include <mpi.h>
#include <cassert>
#include <iostream>
#include <stdint.h>
#include <vector>
#include <omp.h>
#include <string>
#include <sstream>
#include <cstring>
//#include <openacc.h>

#define DENSITY(I,J,K) density_.data_[(I)+(J)*local_nx_+(K)*local_nx_*local_ny_]
#define DENSITY_OLD(I,J,K) density_old_.data_[(I)+(J)*local_nx_+(K)*local_nx_*local_ny_]

//constructor
Diffusion3D::Diffusion3D(val_type dx, count_type nx, count_type ny, count_type nz,coord_type topology,val_type D, val_type T):
						D_(D),T_(T),dx_(dx),
						global_nx_(nx),global_ny_(ny),global_nz_(nz),
						global_lengthx_((nx-1)*dx),global_lengthy_((ny-1)*dx),global_lengthz_((nz-1)*dx)
{
	
	//std::cout << "Constructor Called\n";
	//check # nodes vs mesh points: each node has the same number of meshpoints!
	assert(nx % topology[0] == 0);
	assert(ny % topology[1] == 0);
	assert(nz % topology[2] == 0);
	
	for (unsigned int i =0;i<3;i++)
		topology_[i] = topology[i];
	
	
	dt_ = dx_*dx_*dx_/(6*D);	// dt < dx^2 /(4*D) for stability of FTCS


	
	//determine local # meshpoints, +2 because of ghost cells in each direction
	local_nx_ = global_nx_/topology_[0] + 2;
	local_ny_ = global_ny_/topology_[1] + 2;
	local_nz_ = global_nz_/topology_[2] + 2;
	

	
	setup_MPI_stuff();
	
	density_.resize(local_nx_,local_ny_,local_nz_,0);
	density_old_.resize(local_nx_,local_ny_,local_nz_,0);
	
	set_initial_conditions();
	set_boundary_conditions();
	

}


//set up the communicators, data types etc
void Diffusion3D::setup_MPI_stuff() {
	int periodic[3] = {false,false,false};	//or true,true,true?
	MPI_Cart_create(MPI_COMM_WORLD,3,topology_,periodic,true,&cart_comm_);	//dimensions are the same as in vector
	
	//get rank & size of cartesian communicator
	MPI_Comm_rank(cart_comm_,&rank_);
	MPI_Comm_size(cart_comm_,&size_);
	
	//get coordinates of rank
	MPI_Cart_coords(cart_comm_,rank_,3, cartesian_coords_);
	
	//determine neighbours
	MPI_Cart_shift(cart_comm_,0,1,&left_,&right_);
	MPI_Cart_shift(cart_comm_,1,1,&bottom_,&top_);
	MPI_Cart_shift(cart_comm_,2,1,&back_,&front_);	//z direction points forwards -->back,front (or front,back), not quite sure...
	
	/*std::cout << "Rank: "<< rank_ << " is at position: " << cartesian_coords_[0] << " " << cartesian_coords_[1] << " " << cartesian_coords_[2]
		<< "left right | bottom top | back front " << left_ << " " << right_ << " | " << bottom_ << " " << top_ << " | " << back_ << " " << front_ << "\n";*/
	

	
	//xy-plane is contiguous in memory
	MPI_Type_contiguous(local_nx_*local_ny_,mpi_val_type,&planexy_type_);
	
	//xz-plane: localnz blocks of length localnx with a stride of localnx*localny
	MPI_Type_vector(local_nz_,local_nx_,local_nx_*local_ny_,mpi_val_type,&planexz_type_);
	
	//yz-plane:
/*
	//use indexed type creator
	std::vector<int> offsets(local_ny_*local_nz_);
	int outer_counter = 0;
	#pragma omp parallel for
	for (int k = 0; k < local_nz_;k++,outer_counter++) {
		int inner_counter = 0;
		for (int j = 0; j < local_ny_;j++,inner_counter++) {
			offsets[j+local_ny_*k] = inner_counter*local_nx_ + outer_counter*(local_nx_*local_ny_);
			if (rank_==0)
				std::cout << offsets.at(j+local_ny_*k) << " ";
		}
	} std::cout << "\n";
	std::vector<int> blocklens(local_ny_*local_nz_,1);
	MPI_Type_indexed(local_ny_*local_nz_,blocklens.data(),offsets.data(),mpi_val_type,&planeyz_type_);
	
	Damn I think this is the same as */
	MPI_Type_vector(local_ny_*local_nz_,1,local_nx_,mpi_val_type,&planeyz_type_);

	MPI_Type_commit(&planexy_type_);
	MPI_Type_commit(&planexz_type_);
	MPI_Type_commit(&planeyz_type_);
}


//Destructor
Diffusion3D::~Diffusion3D() {
	//std::cout << "Destructor Called\n";

	//free mpi types
	// leads to errors, as Destructor is called after MPI_Finalize()	
	//MPI_Type_free(&planexy_type_);
	//MPI_Type_free(&planexz_type_);
	//MPI_Type_free(&planeyz_type_);	
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
		<< "\t\t Adress of (0,1,0)\t\t" << &this->density_(0,1,0) <<"\t Stride in y-direction:\t"<<intptr_t(&this->density_(0,1,0))-intptr_t(&this->density_(0,0,0))<<" Bytes\n"
		<< "\t\t Adress of (0,0,1)\t\t" << &this->density_(0,0,1) <<"\t Stride in z-direction:\t"<<intptr_t(&this->density_(0,0,1))-intptr_t(&this->density_(0,0,0))<<" Bytes\n";
	os << "Stride in y-direction should be:\t\t" << density_.get_sizeX()*sizeof(val_type) <<"\n";
	os << "Stried in z-direction should be:\t\t" << density_.get_sizeX()*density_.get_sizeY()*sizeof(val_type)<<"\n";
	os << "\n------------------------------------------------------------------------\n\n";
}

void Diffusion3D::set_boundary_conditions() {}

void Diffusion3D::set_initial_conditions() {
	if (cartesian_coords_[0] == 0 && cartesian_coords_[1] == 0 && cartesian_coords_[2] == 0) {
		#pragma omp parallel for 
		for (int k = 1; k < local_nz_-1;k++) {
			for (int j = 1; j < local_ny_-1;j++) {
				for (int i = 1; i < local_nx_-1;i++) {
					DENSITY(i,j,k) =1;
				}
			}
		}
	}
}

Diffusion3D::val_type Diffusion3D::get_dt() {
	return dt_;
}

void Diffusion3D::start_simulation(count_type stencil) {
	count_type iterations = T_/dt_;
	for (count_type i=0; i<iterations; i++) {
		FTCS();
	}
}

void Diffusion3D::FTCS() {
	// swap densities to densities_old_ vector 
	swap(density_old_,density_);
	
	MPI_Status status[12];
	MPI_Request reqs[12];

        // send left border plane
        MPI_Isend(&density_old_(1,0,0),1,planeyz_type_,left_,0,cart_comm_,&reqs[0]);
        // recieve left ghost plane
        MPI_Irecv(&density_old_(0,0,0),1,planeyz_type_,left_,0,cart_comm_,&reqs[1]);
        // send right border plane
        MPI_Isend(&density_old_(local_nx_-2,0,0),1,planeyz_type_,right_,0,cart_comm_,&reqs[2]);
        // recieve right ghost plane
        MPI_Irecv(&density_old_(local_nx_-1,0,0),1,planeyz_type_,right_,0,cart_comm_,&reqs[3]);
        // send bottom border plane
        MPI_Isend(&density_old_(0,1,0),1,planexz_type_,bottom_,0,cart_comm_,&reqs[4]);
        // recieve bottom ghost plane
        MPI_Irecv(&density_old_(0,0,0),1,planexz_type_,bottom_,0,cart_comm_,&reqs[5]);
        // send top border plane
        MPI_Isend(&density_old_(0,local_ny_-2,0),1,planexz_type_,top_,0,cart_comm_,&reqs[6]);
        // recieve top border plane
        MPI_Irecv(&density_old_(0,local_ny_-1,0),1,planexz_type_,top_,0,cart_comm_,&reqs[7]);
        // send back border plane
        MPI_Isend(&density_old_(0,0,1),1,planexy_type_,back_,0,cart_comm_,&reqs[8]);
        // recieve back ghost plane
        MPI_Irecv(&density_old_(0,0,0),1,planexy_type_,back_,0,cart_comm_,&reqs[9]);
        // send front border plane
        MPI_Isend(&density_old_(0,0,local_nz_-2),1,planexy_type_,front_,0,cart_comm_,&reqs[10]);
        // recieve front border plane
        MPI_Irecv(&density_old_(0,0,local_nz_-1),1,planexy_type_,front_,0,cart_comm_,&reqs[11]);

	// do computation of interior nodes
	val_type prefac = D_*dt_/(dx_*dx_*dx_);
	#pragma omp parallel for schedule(static)
	for (count_type k=2; k<local_nz_-2; k++) {
		for (count_type j=2; j<local_ny_-2; j++) {
			for (count_type i=2;i<local_nx_-2; i++) {
	    			DENSITY(i,j,k) = DENSITY_OLD(i,j,k) + prefac*(	 DENSITY_OLD(i-1,j,k)+DENSITY_OLD(i+1,j,k)
								+DENSITY_OLD(i,j-1,k)+DENSITY_OLD(i,j+1,k)
								+DENSITY_OLD(i,j,k-1)+DENSITY_OLD(i,j,k+1) - 6*DENSITY_OLD(i,j,k) );
			}
		}
	}

	// wait for the ghost cells to arrive
    	MPI_Waitall(12,reqs,status);

	// do computation of border planes
	std::vector<count_type> xLayers, yLayers, zLayers;
	xLayers.push_back(1);
	xLayers.push_back(local_nx_-2);
	yLayers.push_back(1);
	yLayers.push_back(local_ny_-2);
	zLayers.push_back(1);
	zLayers.push_back(local_nz_-2);

	// compute left and right border plane
	#pragma omp parallel for schedule(static)
	for (count_type k=1; k<local_nz_-1; k++) {
		for (count_type j=1; j<local_ny_-1; j++) {
			for (count_type m=0; m<xLayers.size(); m++) {
				count_type i = xLayers[m];
	    			DENSITY(i,j,k) = DENSITY_OLD(i,j,k)+ prefac*(	 DENSITY_OLD(i-1,j,k)+DENSITY_OLD(i+1,j,k)
								+DENSITY_OLD(i,j-1,k)+DENSITY_OLD(i,j+1,k)
								+DENSITY_OLD(i,j,k-1)+DENSITY_OLD(i,j,k+1) - 6*DENSITY_OLD(i,j,k) );
			}		
		}	
	}
	// compute bottom and top border plane
	#pragma omp parallel for schedule(static)
	for (count_type k=1; k<local_nz_-1; k++) {
		for (count_type m=0; m<yLayers.size(); m++) {
			count_type j = yLayers[m];
			for (count_type i=2; i<local_nx_-2; i++) {
	    			DENSITY(i,j,k) = DENSITY_OLD(i,j,k)+ prefac*(	 DENSITY_OLD(i-1,j,k)+DENSITY_OLD(i+1,j,k)
								+DENSITY_OLD(i,j-1,k)+DENSITY_OLD(i,j+1,k)
								+DENSITY_OLD(i,j,k-1)+DENSITY_OLD(i,j,k+1) - 6*DENSITY_OLD(i,j,k) );
			}		
		}	
	}
	// compute left and right border plane
	
	for (count_type m=0; m<zLayers.size(); m++) {
		count_type k = zLayers[m];
		#pragma omp parallel for schedule(static)
		for (count_type j=2; j<local_ny_-2; j++) {
			for (count_type i=2; i<local_nx_-2; i++) {
	    			DENSITY(i,j,k) = DENSITY_OLD(i,j,k)+ prefac*(	 DENSITY_OLD(i-1,j,k)+DENSITY_OLD(i+1,j,k)
								+DENSITY_OLD(i,j-1,k)+DENSITY_OLD(i,j+1,k)
								+DENSITY_OLD(i,j,k-1)+DENSITY_OLD(i,j,k+1) - 6*DENSITY_OLD(i,j,k) );
			}		
		}	
	}
}

void Diffusion3D::print_density(std::string filename) const{
	std::ofstream outfile;
	std::stringstream coord_infox,coord_infoy,coord_infoz,size_infon,size_infom,size_infok;
	coord_infox << cartesian_coords_[0];
	coord_infoy << cartesian_coords_[1];
	coord_infoz << cartesian_coords_[2];
	
	size_infon << local_nx_-2;
	size_infom << local_ny_-2;
	size_infok << local_nz_-2;
	filename = filename + coord_infox.str() +"-" + coord_infoy.str() + "-" +coord_infoz.str() +"_" + size_infon.str() +"-" +size_infom.str() +"-" + size_infok.str() +".txt";
	outfile.open(filename.c_str(),std::ios_base::out);
	for (int k = 1; k < density_.get_sizeZ()-1;k++) {
		for (int j = 1; j < density_.get_sizeY()-1;j++) {
			for (int i = 1; i < density_.get_sizeX()-1;i++) {
				outfile << density_(i,j,k) << " ";
			}
			outfile << "\n";
		}
	}
	
	outfile.close();
	//densities-x-y-z_m_n_k.txt

}


