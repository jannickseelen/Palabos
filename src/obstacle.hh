#ifndef OBSTACLE_HH
#define OBSTACLE_HH

#include "obstacle.h"
#include <palabos3D.hh>
#include "myheaders3D.hh"
#include <string>
#include <exception>

namespace plb{

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Obstacle<T,BoundaryType,SurfaceData,Descriptor>::Obstacle()
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->o.reset(this);
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Obstacle already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Obstacle<T,BoundaryType,SurfaceData,Descriptor>::~Obstacle()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Obstacle";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		objCount--;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Obstacle<T,BoundaryType,SurfaceData,Descriptor>::initialize()
	{
		try{
			std::string meshFileName = Constants<T>::obstacle.fileName;
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Initializing Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] MeshFileName ="+meshFileName;
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			velocity[0] = 0;	velocity[1] = 0;	velocity[2] = 0;
			acceleration[0] = 0;	acceleration[1] = 0;	acceleration[2] = 0;
			rotation[0] = 0;	rotation[1] = 0; rotation[2] = 0;
			rotationalVelocity[0] = 0;	rotationalVelocity[1] = 0;	rotationalVelocity[2] = 0;
			rotationalAcceleration[0] = 0;	rotationalAcceleration[1] = 0;	rotationalAcceleration[2] = 0;
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
			numVertices = triangleSet.getNumVertices();
			numTriangles = triangleSet.getNumTriangles();

			pcout << "The immersed surface  has " << numVertices << " vertices and " << numTriangles << " triangles." << std::endl;
			pcout << "The immersed surface has a maximum triangle edge length of " << maxEdgeLength << std::endl;
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

			for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
				vertices.push_back(triangleSet.getVertex(iVertex));
				T area;
				Array<T,3> unitNormal;
				triangleSet.computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
				areas.push_back(area);
				unitNormals.push_back(unitNormal);
			}

			flowType = voxelFlag::outside;
			volume = getVolume(triangleSet);
			mass = Constants<T>::obstacle.density * volume;
			g = Constants<T>::gravitationalAcceleration;
			velocityFunc.initialize(mass, g, Constants<T>::obstacle.density);

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Volume= "+std::to_string(volume)+" Mass= "+std::to_string(mass);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Done Initializing Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}


	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getCenter(const ConnectedTriangleSet<T>& triangles)
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
			numVertices = triangles.getNumVertices();
			for(int i = 0; i<numVertices; i++){
				Array<T,3> iVertex = triangles.getVertex(i);
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
	Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getCenter()
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
			numVertices = tb->getMesh().getNumVertices();
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
	Box3D Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getDomain(const ConnectedTriangleSet<T>& triangles)
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
	Box3D Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getDomain(const TriangleSet<T>& triangles)
	{
		Box3D d(0,0,0,0,0,0);
		try
		{
			std::vector<Array<Array<T,3>,3> >  iTriangles = triangles.getTriangles();
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

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Box3D Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getDomain()
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
	T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getVolume(const ConnectedTriangleSet<T>& triangles)
	{
		T v = 0;
		try{
			Array<T,3> cg = getCenter(triangles);
			#ifdef PLB_DEBUG
				std::string mesg ="[DEBUG] Calculating Volume";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			numTriangles = triangles.getNumTriangles();
			for(int i =0; i<numTriangles; i++){
				Array<T,3> iTriangle = triangles.getTriangle(i);
				Array<T,3> a = triangleSet.getVertex(iTriangle[0]);
				Array<T,3> b = triangleSet.getVertex(iTriangle[1]);
				Array<T,3> c = triangleSet.getVertex(iTriangle[2]);
				T iVolume = computeTetrahedronSignedVolume(a,b,c,cg);
				if(iVolume<0){iVolume *= -1; }
				v += iVolume;
			}
			if(v<0){
				std::string ex = "[ERROR] Volume= "+std::to_string(v)+" is negative!";
				throw std::runtime_error(ex);
			}
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] DONE Volume= "+std::to_string(v);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return v;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Obstacle<T,BoundaryType,SurfaceData,Descriptor>::moveToStart()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Moving Obstacle to Start Position";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("obstacle").start();
			#endif
				const T dx = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaX();

				Box3D wall_domain = Wall<T,BoundaryType,SurfaceData,Descriptor>::getDomain();
				Array<T,3> wall_cg = Wall<T,BoundaryType,SurfaceData,Descriptor>::getCenter();
				// Find the current location
				Box3D obstacle_domain = getDomain();
				Array<T,3> obstacle_cg = getCenter();

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Obstacle Original Position= "+box_string(obstacle_domain)+" Center= "+array_string(obstacle_cg);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
				T x = 0;
				T y = 0;
				T z = 0;
				x = wall_cg[0] - obstacle_cg[0];
				y = wall_cg[1] - obstacle_cg[1];
				z = wall_domain.z1 - obstacle_domain.z1-1;

				tb->getMesh().translate(Array<T,3>(x,y,z));

				unitNormals.clear();
				unitNormals.resize(numVertices);
				unitNormals.reserve(numVertices);

				areas.clear();
				areas.resize(numVertices);
				areas.reserve(numVertices);

				for(int i = 0; i < numVertices; i++){
					areas[i] = tb->getMesh().computeVertexArea(i);
					unitNormals[i] = tb->getMesh().computeVertexNormal(i);
				}

				obstacle_domain = getDomain();
				obstacle_cg = getCenter();

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Wall Domain= "+box_string(wall_domain)+" Center= "+array_string(wall_cg);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Obstacle Start Position= "+box_string(obstacle_domain)+" Center= "+array_string(obstacle_cg);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] DONE Moving Obstacle to Start Position";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("obstacle").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Obstacle<T,BoundaryType,SurfaceData,Descriptor>::updateImmersedWall()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating Immersed Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("update").start();
			#endif

				numVertices = tb->getMesh().getNumVertices();

				vertices.clear();
				vertices.resize(numVertices);
				vertices.reserve(numVertices);

				unitNormals.clear();
				unitNormals.resize(numVertices);
				unitNormals.reserve(numVertices);

				areas.clear();
				areas.resize(numVertices);
				areas.reserve(numVertices);

				const bool weightedArea = false;

				for(int i = 0; i < numVertices; i++){
					vertices[i] = tb->getMesh().getVertex(i);
					areas[i] = tb->getMesh().computeVertexArea(i);
					unitNormals[i] = tb->getMesh().computeVertexNormal(i,weightedArea);
				}

				//InstantiateImmersedWallData3D<T>(vertices, areas, unitNormals);

				std::vector<MultiBlock3D*> args;
				plint pl = 4;

				args.resize(0);
				args.push_back(Variables<T,BoundaryType,SurfaceData,Descriptor>::container);
				integrateProcessingFunctional(new InstantiateImmersedWallData3D<T>(
					vertices,
					areas,
					unitNormals),
					Variables<T,BoundaryType,SurfaceData,Descriptor>::container->getBoundingBox(),
					*Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice, args, pl);


			#ifdef PLB_DEBUG
				mesg =   "[DEBUG] DONE Updating Immersed Wall";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("update").stop();
			#endif
			}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	bool Obstacle<T,BoundaryType,SurfaceData,Descriptor>::move()
	{
		bool stop = false;
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Moving Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("move").start();
			#endif
				const T dt = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaT();
				const T dx = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaX();
				const T omega = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getOmega();
				const T rho_LB = (T)1.0;
				const T timeLB = Variables<T,BoundaryType,SurfaceData,Descriptor>::time;
				const Box3D lattice_domain = Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBoundingBox();
				const Box3D obstacle_domain = getDomain();
				normalFunc.update(tb.get());

				T factor = util::sqr(util::sqr(dx)) / util::sqr(dt);

				for (int i = 0; i < Constants<T>::ibIter; i++){
					resetForceStatistics<T>(*Variables<T,BoundaryType,SurfaceData,Descriptor>::container);
				}

				recomputeImmersedForce<T>(normalFunc, omega, rho_LB,
					*Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice,
					*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
					Constants<T>::envelopeWidth, obstacle_domain, true);

				Array<T,3> force = Array<T,3>(0,0,0);
				force = -reduceImmersedForce<T>(*Variables<T,BoundaryType,SurfaceData,Descriptor>::container);

				Array<T,3> center = getCenter();

				Array<T,3> torque = Array<T,3>(0,0,0);
				torque = -reduceAxialTorqueImmersed(*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
										center, Array<T,3>(1,1,1));

				stop = velocityFunc.update(Variables<T,BoundaryType,SurfaceData,Descriptor>::p,
											timeLB,force,torque,tb.get(),lattice_domain);
				/*
				for (int i = 0; i < Constants<T>::ibIter; i++){
					indexedInamuroIteration<T>(velocityFunc,
									*Variables<T,BoundaryType,SurfaceData,Descriptor>::rhoBar,
									*Variables<T,BoundaryType,SurfaceData,Descriptor>::j,
									*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
									Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getTau(),
									true);
				}*/

				normalFunc.update(tb.get());

				updateImmersedWall();

			#ifdef PLB_DEBUG
				mesg =   "[DEBUG] DONE Moving Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("move").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return stop;
	}


} // namespace plb

#endif //OBSTACLE_HH

