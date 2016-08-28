#ifndef OBSTACLE_HH
#define OBSTACLE_HH

#include "obstacle.h"
#include <Palabos3D.hh>
#include "myheaders3D.hh"
#include <string>
#include <exception>

namespace plb{

	template<typename T, class BoundaryType, template<class U> class Descriptor>
	Obstacle<T, BoundaryType, Descriptor>::Obstacle()
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

	template<typename T, class BoundaryType, template<class U> class Descriptor>
	Obstacle<T, BoundaryType, Descriptor>::~Obstacle()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Obstacle";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		objCount--;
	}

	template<typename T, class BoundaryType, template<class U> class Descriptor>
	void Obstacle<T, BoundaryType, Descriptor>::initialize()
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
			T x = Constants<T>::obstacle_data[0];
			T y = Constants<T>::obstacle_data[1];
			T z = Constants<T>::obstacle_data[2];
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
					triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
					global::mpiData().sendTriangleSet<T>(triangleSet);
				}
				else{ triangleSet = global::mpiData().receiveTriangleSet<T>(); }
			#else
				triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
			#endif
			flowType = voxelFlag::outside;
			volume = getVolume();
			mass = density * volume;
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Number of triangles in Mesh = "+std::to_string(triangleSet.getTriangles().size());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Done Initializing Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){
			int line = __LINE__;
			std::string file = __FILE__;
			std::string what = e.what();
			std::string ex = "[ERROR]: "+what+" [FILE:"+file+",LINE:"+std::to_string(line)+"]";
			global::log(ex);
			std::cerr << ex <<std::endl;
			throw e;
		}
	}


	template<typename T, class BoundaryType, template<class U> class Descriptor>
	Obstacle<T, BoundaryType, Descriptor>& Obstacle<T, BoundaryType, Descriptor>::getCenter()
	{
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

	template<typename T, class BoundaryType, template<class U> class Descriptor>
	T Obstacle<T, BoundaryType, Descriptor>::getVolume(){
		T volume = 0;
		getCenter();
		std::vector<Array<Array<T,3>,3> > triangles = triangleSet.getTriangles();
		for(int i =0; i<triangles.size(); i++){
			Pyramid<T> p(triangles[i],center);
			volume += p.volume();
		}
		return volume;
	}

	template<typename T, class BoundaryType, template<class U> class Descriptor>
	void Obstacle<T, BoundaryType, Descriptor>::move(const plint& dt, std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice)
	{
		Array<T,3> coord = mesh->getPhysicalLocation();
		Array<T,3> fluidForce = bc->getForceOnObject();
		acceleration[0] = fluidForce[0]/mass;
		acceleration[1] = fluidForce[1]/mass;
		acceleration[2] = fluidForce[2]/mass;
		acceleration[2] -= Constants<T>::gravitationalAcceleration;
		velocity[0] += acceleration[0]*dt;
		velocity[1] += acceleration[1]*dt;
		velocity[2] += acceleration[2]*dt;
		Array<T,3> ds(velocity[0]*dt, velocity[1]*dt, velocity[2]*dt);
		coord[0] += ds[0]; coord[1] += ds[1]; coord[2] += ds[2];
		mesh->setPhysicalLocation(coord);
		reparallelize(*lattice);
	}

	template<typename T, class BoundaryType, template<class U> class Descriptor>
	void Obstacle<T, BoundaryType, Descriptor>::move()
	{
		Array<T,3> vec;
		vec[0] = 0;
		vec[1] = 0;
		vec[2] = 0;
		triangleSet.translate(vec);
	}

} // namespace plb

#endif //OBSTACLE_HH

