//HPCSE PROJECT
//Stefan H & Maurice G
//3d Diffusion solver

#ifndef DIFF3D_HPP 
#define DIFF3D_HPP 

#include <mpi.h>
#include <vector>
#include <cassert>
#include <string>
#include <iostream>
#include <sstream>
#include <fstream>


/*
* \class Diffusion3D
* \brief Class implementing a (distributed,mesh-based) 3d Diffusion solver on a box shaped domain
*/

class Diffusion3D {
 
 	public:
 	

		typedef double val_type;
		typedef std::size_t count_type;
		typedef std::vector<val_type> domain3d_type; 		//todo: create a better type using aligned allocators, s.t. A(1,2,3) is possible etc.
		typedef int coord_type[3];						//std::array is c++11 --> no cray/pgi compiler support
		#define mpi_val_type MPI_DOUBLE
	
		//todo: public functions: Constructor, Destructor, Run, write..
		
	 	/*!
	 	* Constructor,
	 	* @param dx is the meshwidth of the mesh (uniform in all directions)
	 	* @param nx is the number of meshpoints in x direction
	 	* @param ny is the number of meshpoints in y direction 
	 	* @param nz is the number of meshpoints in z direction 
	 	* @param a is a corner of the rectangular domain. a&b together define the domain
	 	* @param b is a corner of the rectangular domain
	 	* @param D is the diffusion constant to be used
	 	* @param T is the wanted end-time of the simulation
	 	*/
		Diffusion3D(val_type dx, val_type nx, val_type ny, val_type nz,
					coord_type a, coord_type b,
					val_type D, val_type T);
					
		/*!
		* Destructor, frees mpi datatypes etc
		*/
		~Diffusion3D();
		
		/*!
		* This function starts the simulation and runs it until time T
		* @param stencil specifies what stencil(FTCS, BTCS,...) to use: 1 stands for..., 2 for...
		*/
		void start_simulation(count_type stencil);
		
		/*!
		* This function writes the content of density into a file
		* @param filename is the name of the desired output-file
		*/
		void print_density (std::string filename);
		
		
		
		
		
		
		
	private:
		//--------------------everything here is the same for each rank---------------//
		//problem specific parameters
		val_type 	D_;					// Diffusion coefficient
		val_type 	T_;					// Time to simulate
	
		
		
		// numerical parameters
		val_type 	dt_; 				// timestep
		val_type 	dx_;				// meshwidth, uniform in all directions
		
		
		
		//(global) domain specific parameters
		val_type 	global_lengthx_,
				 	global_lengthy_,
				 	global_lengthz_; 	// size of global domain
		
		count_type  global_nx_,
					global_ny_,
					global_nz_;  		//# mesh points in each direction (global)
					
		
		
		//---------------------------variables which have rank specific  values -------------
		
		domain3d_type density,			//the densities of each gridpoint are stored in here
					  density_old;
		//MPI Stuff
		int 		rank_;
		int 		size_;
		int 		cartesian_coords[3];	
		
		count_type	local_nx_,
					local_ny_,
					local_nz_;			// #mesh points in each direction (local)
					
		int 		left_,
					right_,
					top_,
					bottom_,
					front_,
					back_;				// neighbours in cartesian mesh 
					
		MPI_Comm 	cart_comm_;			//cartesian communicator
		
		MPI_Datatype planexy_type_,		//datatypes for ghost planes
					 planexz_type_,
					 planeyz_type_;		
		
		//maybe add MPI_Request, MPI_Status -->don't create them every iteration
		
		//-----------------------private functions------------------------
		//todo: private functions like initial config, b.c. etc








}


#endif