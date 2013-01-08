//HPCSE PROJECT
//Stefan H & Maurice G
//3d Diffusion solver

#ifndef DIFF3D_HPP 
#define DIFF3D_HPP 

#include <mpi.h>
#include <vector>

/*
* \class Diffusion3D
* \brief Class implementing a (distributed,mesh-based) 3d Diffusion solver on a box shaped domain
*/

class Diffusion3D {
 
 	public:
 	

		typedef double val_type;
		typedef std::size_t count_type;
		typedef std::vector<val_type> domain_type; 		//todo: create a better type using aligned allocators, s.t. A(1,2,3) is possible etc.
		typedef std::array<val_type,3> coord_type;
		#define mpi_val_type MPI_DOUBLE
	
		//todo: public functions: Constructor, Destructor, Run, write..
		
		
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
					
		
		
		//---------------------------variables have seperate values on each rank--------------
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
					
		//stuff todo: cartesian communicator, types for the (ghost)planes, status,requests
					
		
		
		//-----------------------private functions------------------------
		//todo: private functions like initial config, b.c. etc








}


#endif
