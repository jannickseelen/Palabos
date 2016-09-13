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
			std::string meshFileName = Constants<T>::c->obstacle_file;
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Initializing Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] MeshFileName ="+meshFileName;
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			dynamicMesh = Constants<T>::dynamicObstacle;
			T x = Constants<T>::obstacle_data[0];
			T y = Constants<T>::obstacle_data[1];
			T z = Constants<T>::obstacle_data[2];
			location = Array<T,3>(x,y,z);
			position = Point<T>(x,y,z);
			referenceDirection = Constants<T>::obstacle_data[3];
			density = Constants<T>::obstacle_data[4];

			velocity[0] = 0;	velocity[1] = 0;	velocity[2] = 0;
			acceleration[0] = 0;	acceleration[1] = 0;	acceleration[2] = 0;
			rotation[0] = 0;	rotation[1] = 0; rotation[2] = 0;
			rotationalVelocity[0] = 0;	rotationalVelocity[1] = 0;	rotationalVelocity[2] = 0;
			rotationalAcceleration[0] = 0;	rotationalAcceleration[1] = 0;	rotationalAcceleration[2] = 0;
			#ifdef PLB_MPI_PARALLEL
				if(global::mpi().isMainProcessor()){
					triangleSet = TriangleSet<T>(meshFileName, Constants<T>::precision, STL);
					global::mpiData().sendTriangleSet<T>(triangleSet);
				}
				else{ triangleSet = global::mpiData().receiveTriangleSet<T>(); }
			#else
				triangleSet = TriangleSet<T>(meshFileName, Constants<T>::precision, STL);
			#endif
			flowType = voxelFlag::outside;
			getVolume();
			pcout << "VOLUME= " << std::to_string(volume) << std::endl;
			mass = density * volume;
			T g = Constants<T>::gravitationalAcceleration;
			surfaceVelocity.initialize(location, mass, g);
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Number of triangles in Mesh = "+std::to_string(triangleSet.getTriangles().size());
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
	Obstacle<T,BoundaryType,SurfaceData,Descriptor>& Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getCenter()
	{
		try{
			Cuboid<T> cuboid = triangleSet.getBoundingCuboid();
			Array<T,3>	lowerLeftCorner = cuboid.lowerLeftCorner;
			Array<T,3>	upperRightCorner = cuboid.upperRightCorner;
			T lowerBound = 0;
			T upperBound = 0;
			lowerBound = std::min(lowerLeftCorner[0],upperRightCorner[0]);
			upperBound = std::max(lowerLeftCorner[0],upperRightCorner[0]);
			center.x = (upperBound-lowerBound)/2;
			lowerBound = std::min(lowerLeftCorner[1],upperRightCorner[1]);
			upperBound = std::max(lowerLeftCorner[1],upperRightCorner[1]);
			center.y = (upperBound-lowerBound)/2;
			lowerBound = std::min(lowerLeftCorner[2],upperRightCorner[2]);
			upperBound = std::max(lowerLeftCorner[2],upperRightCorner[2]);
			center.z = (upperBound-lowerBound)/2;
			return *o;
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getVolume(){
		try{
			getCenter();
			std::vector<Array<Array<T,3>,3> > triangles = triangleSet.getTriangles();
			for(int i =0; i<triangles.size(); i++){
				Pyramid<T> p(triangles[i],center);
				volume += p.volume();
			}
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
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
				mesh->getMesh().translate(location);
				force.reset(new GetForceOnObjectFunctional3D<T,BoundaryType>(model->clone()));
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] DONE Moving Obstacle to Start Position";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("obstacle").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getForce()
	{
		Array<T,3> f = Array<T,3>(0,0,0);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Calculating Force on Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("obstacle").start();
			#endif
				force.reset(new GetForceOnObjectFunctional3D<T,BoundaryType>(model->clone()));

				PlainReductiveBoxProcessingFunctional3D* func = nullptr;
				func = dynamic_cast<PlainReductiveBoxProcessingFunctional3D*>(force.get());

				Box3D domain = bc->getArg().getBoundingBox();

				std::map<plint,BlockLattice3D<T,Descriptor>*> blocks =
					Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBlockLattices();
				std::vector<AtomicBlock3D*> atomics;
				plint size = blocks.size();
				atomics.resize(size);
				atomics.reserve(size);
				for(int i = 0; i<size; i++){
					atomics.push_back(dynamic_cast<AtomicBlock3D*>(blocks[i]));
				}

				if(force != nullptr){applyProcessingFunctional(*func, domain, atomics); }
				else{ throw std::runtime_error("GetForceOnObjectFunctional3D not Properly Initialized"); }

				f += force->getForce();
			#ifdef PLB_DEBUG
				mesg =  "[DEBUG] DONE Calculating Force on Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("obstacle").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return f;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Obstacle<T,BoundaryType,SurfaceData,Descriptor>::move()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Moving Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("move").start();
			#endif
				const T dt = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaT();
				const T dx = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaX();
				Array<T,3> force = Array<T,3>(0,0,0);
				force = getForce();
				std::vector<Array<T,3> > vertexList = mesh->getVertexList();
				Array<T,3> ds = Array<T,3>(0,0,0);
				ds = surfaceVelocity.update(Variables<T,BoundaryType,SurfaceData,Descriptor>::time,force,dt,dx);
				if(ds[0] != 0 && ds[1] != 0 && ds[2] != 0){
					for(int i = 0; i<vertexList.size(); i++){
						vertexList[i] += ds;
					}
					instantiateImmersedWallData(mesh->getVertexList(),
												mesh->getAreaList(),
												*Variables<T,BoundaryType,SurfaceData,Descriptor>::container);
					for (int i = 0; i < Constants<T>::ibIter; i++){
						inamuroIteration<T>(*velocityFunc,
										*Variables<T,BoundaryType,SurfaceData,Descriptor>::rhoBar,
										*Variables<T,BoundaryType,SurfaceData,Descriptor>::j,
										*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
										Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getTau(),
										true);
					}
				}
				else{ throw std::runtime_error("Could not move obstacle since its new position is the same as the old one"); }
			#ifdef PLB_DEBUG
				mesg =   "[DEBUG] Moving Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("move").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	

} // namespace plb

#endif //OBSTACLE_HH

