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
			std::string meshFileName = Constants<T>::wall.fileName;
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

			Box3D domain = getDomain(surface);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Domain BEFORE Scaling "+ box_string(domain) +" in physical units";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			T x = 0;
			T y = 0;
			T z = 0;
			if(domain.x0<0){ x = -domain.x0;}
			if(domain.y0<0){ y = -domain.y0;}
			if(domain.z0<0){ z = -domain.z0;}
			Array<T,3> shift = Array<T,3>(x,y,z);
			surface.translate(shift);

			T alpha = 0.01;
			surface.scale(alpha);
			
			domain = getDomain(surface);

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Domain AFTER scaling "+ box_string(domain) +" in physical units";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			const T dx = Constants<T>::lb.dx;

			T maxEdgeLength = surface.getMaxEdgeLength();

			triangleSet = ConnectedTriangleSet<T>(surface);
			plint numVertices = triangleSet.getNumVertices();
			plint numTriangles = triangleSet.getNumTriangles();

			pcout << "The wall surface  has " << numVertices << " vertices and " << numTriangles << " triangles." << std::endl;
			pcout << "The wall surface has a maximum triangle edge length of " << maxEdgeLength << std::endl;
			pcout << "dx = "<<dx <<std::endl;

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
			domain = getDomain(triangleSet);
			#ifdef PLB_DEBUG
				mesg="[DEBUG] Done Initializing Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Array<T,3> Wall<T,BoundaryType,SurfaceData,Descriptor>::getCenter()
	{
		Array<T,3> cg = Array<T,3>(0,0,0);
		try{
			#ifdef PLB_DEBUG
				std::string mesg ="[DEBUG] Calculating Center";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			T x = 0;
			T y = 0;
			T z = 0;
			T numVertices = tb->getMesh().getNumVertices();
			for(int i = 0; i<numVertices; i++){
				Array<T,3> iVertex = tb->getMesh().getVertex(i);
				x += iVertex[0];
				y += iVertex[1];
				z += iVertex[2];
			}
			cg = Array<T,3>(x/numVertices, y/numVertices, z/numVertices);
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] DONE Center= "+array_string(cg);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return cg;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Box3D Wall<T,BoundaryType,SurfaceData,Descriptor>::getDomain()
	{
		Box3D d(0,0,0,0,0,0);
		try
		{
			T numVertices = tb->getMesh().getNumVertices();
			T zmin = std::numeric_limits<T>::max();
			T zmax = std::numeric_limits<T>::min();
			T ymin = std::numeric_limits<T>::max();
			T ymax = std::numeric_limits<T>::min();
			T xmin = std::numeric_limits<T>::max();
			T xmax = std::numeric_limits<T>::min();
			for(int i = 0; i<numVertices; i++){
				Array<T,3> iVertex = tb->getMesh().getVertex(i);
				if(iVertex[0] < xmin){ xmin = iVertex[0]; }
				if(iVertex[0] > xmax){ xmax = iVertex[0]; }
				if(iVertex[1] < ymin){ ymin = iVertex[1]; }
				if(iVertex[1] > ymax){ ymax = iVertex[1]; }
				if(iVertex[2] < zmin){ zmin = iVertex[2]; }
				if(iVertex[2] > zmax){ zmax = iVertex[2]; }
			}
			d = Box3D(xmin,xmax,ymin,ymax,zmin,zmax);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return d;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Box3D Wall<T,BoundaryType,SurfaceData,Descriptor>::getDomain(const ConnectedTriangleSet<T>& triangles)
	{
		Box3D d(0,0,0,0,0,0);
		try
		{
			T numVertices = triangles.getNumVertices();
			T zmin = std::numeric_limits<T>::max();
			T zmax = std::numeric_limits<T>::min();
			T ymin = std::numeric_limits<T>::max();
			T ymax = std::numeric_limits<T>::min();
			T xmin = std::numeric_limits<T>::max();
			T xmax = std::numeric_limits<T>::min();
			for(int i = 0; i<numVertices; i++){
				Array<T,3> iVertex = triangles.getVertex(i);
				if(iVertex[0] < xmin){ xmin = iVertex[0]; }
				if(iVertex[0] > xmax){ xmax = iVertex[0]; }
				if(iVertex[1] < ymin){ ymin = iVertex[1]; }
				if(iVertex[1] > ymax){ ymax = iVertex[1]; }
				if(iVertex[2] < zmin){ zmin = iVertex[2]; }
				if(iVertex[2] > zmax){ zmax = iVertex[2]; }
			}
			d = Box3D(xmin,xmax,ymin,ymax,zmin,zmax);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return d;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Box3D Wall<T,BoundaryType,SurfaceData,Descriptor>::getDomain(const TriangleSet<T>& triangles)
	{
		Box3D d(0,0,0,0,0,0);
		try
		{
			std::vector<Array<Array<T,3>,3> > iTriangles = triangles.getTriangles();
			T zmin = std::numeric_limits<T>::max();
			T zmax = std::numeric_limits<T>::min();
			T ymin = std::numeric_limits<T>::max();
			T ymax = std::numeric_limits<T>::min();
			T xmin = std::numeric_limits<T>::max();
			T xmax = std::numeric_limits<T>::min();
			for(int i = 0; i<iTriangles.size(); i++){
				Array<Array<T,3>,3> iTriangle = iTriangles[i];
				for(int v = 0; v<3; v++){
					Array<T,3> iVertex = iTriangle[v];
					if(iVertex[0] < xmin){ xmin = iVertex[0]; }
					if(iVertex[0] > xmax){ xmax = iVertex[0]; }
					if(iVertex[1] < ymin){ ymin = iVertex[1]; }
					if(iVertex[1] > ymax){ ymax = iVertex[1]; }
					if(iVertex[2] < zmin){ zmin = iVertex[2]; }
					if(iVertex[2] > zmax){ zmax = iVertex[2]; }
				}
			}
			d = Box3D(xmin,xmax,ymin,ymax,zmin,zmax);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return d;
	}

}// namespace plb

#endif //WALL_H

