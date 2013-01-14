//Stefan H
//(Heavily influenced by "matrix.hpp" hpc12, Andreas Hehn)

#ifndef VEC3_HPP
#define VEC3_HPP

#include <vector>
#include "aligned_allocator.hpp"
#include <cassert>

namespace D3 {

/*
* \class Vector3d
* \brief Class implementing a 3d-vector (possibly aligned), which grants acces with (,,) operator
*/
template <typename T,typename Allocator = hpc12::aligned_allocator<T,64> >
class Vector3d {
	public:
		typedef T val_type;
		
		/*!
		* Constructor
		* @param x is the size of the domain in x-direction 
		* @param y is the size of the domain in y-direction 
		* @param z is the size of the domain in y-direction 
		*/
		explicit Vector3d (unsigned int x=0, unsigned int y=0, unsigned int z=0) //default parameters (0) are necessary,else compiler complains when calling the
				 : data_(x*y*z,0),												 //constructor of a class where Vector3d is a member
				 size_X_(x),size_Y_(y),size_Z_(z)
				 {}
		
		/*!
		* Element Access: Elements are ordered : continuous in x-direction (no stride), 
												 then y-direction(stride: size_X_),
												 then z-direction(stride: size_X_*size_Y_)
												 
		* @param i -> x-direction 
		* @param j -> y-direction 
		* @param k -> z-direction 
		*/
		val_type& operator()(unsigned int i, unsigned int j, unsigned int k) {
					assert(i < size_X_);
					assert(j < size_Y_);
					assert(k < size_Z_);
					return data_[i + size_X_*j + size_X_*size_Y_*k];
					}
		
		val_type const& operator()(unsigned int i, unsigned int j, unsigned int k) const {
					assert(i < size_X_);
					assert(j < size_Y_);
					assert(k < size_Z_);
					return data_[i + size_X_*j + size_X_*size_Y_*k];
					}
		
		unsigned int get_sizeX () const { return size_X_;} 
		unsigned int get_sizeY () const { return size_Y_;}
		unsigned int get_sizeZ () const { return size_Z_;}
				 
		void resize(unsigned int Nx,unsigned int Ny, unsigned int Nz, val_type value=0) { 
				assert(Nx >0);
				assert(Ny > 0);
				assert(Nz > 0);
				data_.resize(Nx*Ny*Nz,value);
				size_X_ = Nx;
				size_Y_ = Ny;
				size_Z_ = Nz;
		}
		
		friend void swap(Vector3d & A, Vector3d & B) {
			using std::swap;
			swap(A.data_,B.data_);
			swap(A.size_X_,B.size_X_);
			swap(A.size_Y_,B.size_Y_);
			swap(A.size_Z_,B.size_Z_);
		}
				 
	
	private:
		std::vector<val_type,Allocator> data_;
		unsigned int size_X_,size_Y_,size_Z_;
		
};

//  Output x,y-plane(k=0) first, than k++
template <class T, class Allocator>
std::ostream & operator << (std::ostream &os, Vector3d<T,Allocator>  const& v) {
	os <<"#           		^                    \n" // after trial and error...
	   <<"#		    (y) |	                 \n"
	   <<"#			|  /				 \n"
	   <<"#			| / 				 \n"
	   <<"#			|/					 \n"
	   <<"#		  ------/----------> (x)     \n"
	   <<"#		       /|					 \n"
	   <<"#                     / |					 \n"
	   <<"#                    /  |					 \n"
	   <<"#                   /   |					 \n"
	   <<"#	    (z)    v                    \n\n\n";
		
	for (unsigned int k = 0; k < v.get_sizeZ();k++) {
		os << "Next level: k=" << k <<"\n";
		for (unsigned int j = 0; j < v.get_sizeY();j++) {
			for (unsigned int i = 0; i < v.get_sizeX();i++) {
				if (i == 0) os << "\n";
				if (i != 0) os << " ";
				os << v(i,j,k);
			}
			os << "\n";
		}
		os << "\n";
	}
	return os;
}


}
#endif
