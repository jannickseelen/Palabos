#ifndef WALL_HH
#define WALL_HH

#include "wall.h"
#include <palabos3D.hh>
#include "myheaders3D.hh"

namespace plb{

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Wall<T,BoundaryType,SurfaceData,Descriptor>::Wall()
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->w.reset(this);
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Wall already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Wall<T,BoundaryType,SurfaceData,Descriptor>::~Wall()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Wall";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		objCount--;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Wall<T,BoundaryType,SurfaceData,Descriptor>::initialize()
	{
		try
		{
			std::string meshFileName = Constants<T>::c->wall_file;
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] MeshFileName ="+meshFileName;
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			#ifdef PLB_MPI_PARALLEL
				if(global::mpi().isMainProcessor()){
					triangleSet = TriangleSet<T>(meshFileName, LDBL, STL);
					global::mpiData().sendTriangleSet<T>(triangleSet);}
				else{ triangleSet = global::mpiData().receiveTriangleSet<T>(); }
			#else
				triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
			#endif
			flowType = voxelFlag::inside;
			temperature = Constants<T>::wall_data[1];
			referenceDirection = Constants<T>::wall_data[0];
			#ifdef PLB_DEBUG
				mesg="[DEBUG] Number of triangles in Mesh = "+std::to_string(triangleSet.getTriangles().size());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg="[DEBUG] Done Initializing Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}


}// namespace plb

#endif //WALL_H

