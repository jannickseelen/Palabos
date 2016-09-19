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
			TriangleSet<T> surface;
			#ifdef PLB_MPI_PARALLEL
				if(global::mpi().isMainProcessor()){
					surface = TriangleSet<T>(meshFileName, Constants<T>::precision, STL);
					global::mpiData().sendTriangleSet<T>(surface);
				}
				else{ surface = global::mpiData().receiveTriangleSet<T>(); }
			#else
				surface = TriangleSet<T>(meshFileName, Constants<T>::precision, STL);
			#endif

			const T dx = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaX();
			T maxEdgeLength = surface.getMaxEdgeLength();
			surface.scale(1.0 / dx);
			surface.translate(	location / dx);

			triangleSet = ConnectedTriangleSet<T>(surface);
			plint numVertices = triangleSet.getNumVertices();
			plint numTriangles = triangleSet.getNumTriangles();

			pcout << "The wall surface  has " << numVertices << " vertices and " << numTriangles << " triangles." << std::endl;
			pcout << "The wall surface has a maximum triangle edge length of " << maxEdgeLength << std::endl;

			if (maxEdgeLength >= 4.0 * dx) {
				pcout << std::endl;
				pcout << "CAUTION: The maximum triangle edge length for the immersed surface is greater than "
					  << " 4 times dx."
					  << std::endl;
				pcout << "         The immersed boundary method will not work correctly. Surface refinement is necessary."
					  << std::endl;
				pcout << std::endl;
				exit(1);
			} else if (maxEdgeLength > dx) {
				pcout << std::endl;
				pcout << "WARNING: The maximum triangle edge length for the immersed surface is greater than dx."
					  << std::endl;
				pcout << "         The immersed boundary method might not work in an optimal way. Surface refinement is recommended."
					  << std::endl;
				pcout << std::endl;
			}
			flowType = voxelFlag::inside;
			temperature = Constants<T>::wall_data[1];
			referenceDirection = Constants<T>::wall_data[0];
			dynamicMesh = Constants<T>::dynamicWall;
			#ifdef PLB_DEBUG
				mesg="[DEBUG] Done Initializing Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}


}// namespace plb

#endif //WALL_H

