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

			Cuboid<T> cube = surface.getBoundingCuboid();
			Array<T,3> lowerLeftCorner = cube.lowerLeftCorner;
			Array<T,3> upperRightCorner = cube.upperRightCorner;

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Bounded Cuboid BEFORE Scaling Lower Left Corner "+array_string(lowerLeftCorner)+" Upper Right Corner "+
					array_string(upperRightCorner)+" in physical units";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			T x = 0;
			T y = 0;
			T z = 0;
			if(lowerLeftCorner[0]<0){ x = -lowerLeftCorner[0];}
			if(lowerLeftCorner[1]<0){ y = -lowerLeftCorner[1];}
			if(lowerLeftCorner[2]<0){ z = -lowerLeftCorner[2];}
			Array<T,3> shift = Array<T,3>(x,y,z);
			surface.translate(shift);

			//T xLength = upperRightCorner[0] - lowerLeftCorner[0];
			//T alpha = Constants<T>::obstacle.dim[0] / xLength;
			T alpha = 0.01;
			surface.scale(alpha);

			cube = surface.getBoundingCuboid();
			lowerLeftCorner = cube.lowerLeftCorner;
			upperRightCorner = cube.upperRightCorner;

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Bounded Cuboid AFTER Scaling Lower Left Corner "+array_string(lowerLeftCorner)+" Upper Right Corner "+
					array_string(upperRightCorner)+" in physical units";
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
			getVolume();
			pcout << "VOLUME= " << std::to_string(volume) << std::endl;
			mass = Constants<T>::obstacle.density * volume;
			g = Constants<T>::gravitationalAcceleration;

			if(SurfaceVelocity<T>::objCount==0){ velocityFunc = SurfaceVelocity<T>(); }
			if(SurfaceNormal<T>::objCount==0){ normalFunc = SurfaceNormal<T>(); }
			normalFunc.update(triangleSet);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Done Initializing Obstacle";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}


	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getCenter()
	{
		Array<T,3> cg = Array<T,3>(0,0,0);
		try{
			T x = 0;
			T y = 0;
			T z = 0;
			numVertices = triangleSet.getNumVertices();
			for(int i = 0; i<numVertices; i++){
				Array<T,3> iVertex = vertices[i];
				x += iVertex[0];
				y += iVertex[1];
				z += iVertex[2];
			}
			cg = Array<T,3>(x/numVertices, y/numVertices, z/numVertices);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return cg;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getVolume(){
		T v = 0;
		try{
			Array<T,3> cg = getCenter();
			numTriangles = triangleSet.getNumTriangles();
			for(int i =0; i<numTriangles; i++){
				Array<T,3> iTriangle = triangleSet.getTriangle(i);
				Array<T,3> a = triangleSet.getVertex(iTriangle[0]);
				Array<T,3> b = triangleSet.getVertex(iTriangle[1]);
				Array<T,3> c = triangleSet.getVertex(iTriangle[2]);
				T iVolume = computeTetrahedronSignedVolume(a,b,c,cg);
				v += iVolume;
			}
			volume = v;
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
				Array<T,3> wall_cg = Wall<T,BoundaryType,SurfaceData,Descriptor>::center;
				// Find the current location
				T x = 0;
				T y = 0;
				T z = 0;
				T zmin = 0;
				T zmax = 0;
				T ymin = 0;
				T ymax = 0;
				T xmin = 0;
				T xmax = 0;
				numVertices = triangleSet.getNumVertices();
				for(int i = 0; i<numVertices; i++){
					Array<T,3> iVertex = triangleSet.getVertex(i);
					x += iVertex[0];
					if(iVertex[0] < xmin){ xmin = iVertex[0]; }
					if(iVertex[0] > xmax){ xmax = iVertex[0]; }
					y += iVertex[1];
					if(iVertex[1] < ymin){ ymin = iVertex[1]; }
					if(iVertex[1] > ymax){ ymax = iVertex[1]; }
					z += iVertex[2];
					if(iVertex[2] < zmin){ zmin = iVertex[2]; }
					if(iVertex[2] > zmax){ zmax = iVertex[2]; }
				}
				Array<T,3> center = Array<T,3>(x/numVertices, y/numVertices, z/numVertices);
				Box3D obstacle_domain = Box3D(xmin,xmax,ymin,ymax,zmin,zmax);

				x = 0;
				y = 0;
				z = 0;
				x = wall_cg[0] - center[0];
				y = wall_cg[1] - center[1];
				z = wall_domain.z1 - zmax;
				TriangleSet<T> simple = *triangleSet.toTriangleSet(Constants<T>::precision);
				simple.translate(Array<T,3>(x,y,z));
				triangleSet = ConnectedTriangleSet<T>(simple);

				// Move to new location
				//std::vector<Array<T,3> > newVertices;
				//newVertices.resize(numVertices);
				//newVertices.reserve(numVertices);

				unitNormals.clear();
				unitNormals.resize(numVertices);
				unitNormals.reserve(numVertices);

				areas.clear();
				areas.resize(numVertices);
				areas.reserve(numVertices);

				for(int i = 0; i < numVertices; i++){
					/*Array<T,3> iVertex = triangleSet.getVertex(i);
					iVertex[0] += x;
					if(iVertex[0] < xmin){ xmin = iVertex[0]; }
					if(iVertex[0] > xmax){ xmax = iVertex[0]; }
					iVertex[1] += y;
					if(iVertex[1] < ymin){ ymin = iVertex[1]; }
					if(iVertex[1] > ymax){ ymax = iVertex[1]; }
					iVertex[2] += z;
					if(iVertex[2] < zmin){ zmin = iVertex[2]; }
					if(iVertex[2] > zmax){ zmax = iVertex[2]; }
					newVertices[i] = iVertex;*/
					T area = 0;
					Array<T,3> unitNormal = Array<T,3>(0,0,0);
					triangleSet.computeVertexAreaAndUnitNormal(i, area, unitNormal);
					unitNormals[i] = unitNormal;
					areas[i] = area;
				}
				//triangleSet.swapGeometry(newVertices);

				Box3D lattice = Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBoundingBox();

				if(lattice.x0 > obstacle_domain.x0 || lattice.x1 < obstacle_domain.x1
				|| lattice.y0 > obstacle_domain.y0 || lattice.y1 < obstacle_domain.y1
				|| lattice.z0 > obstacle_domain.z0 || lattice.z1 < obstacle_domain.z1
				|| lattice.x0 > wall_domain.x0 || lattice.x1 < wall_domain.x1
				|| lattice.y0 > wall_domain.y0 || lattice.y1 < wall_domain.y1
				|| lattice.z0 > wall_domain.z0 || lattice.z1 < wall_domain.z1)
				{
					std::string ex = "[ERROR] Domain mismatch Lattice = "+ box_string(lattice)+" but obstacle = "
					+box_string(obstacle_domain)+" and wall = "+box_string(wall_domain);
					throw std::runtime_error(ex);
				}

				instantiateImmersedWallData(vertices, areas, unitNormals,	*Variables<T,BoundaryType,SurfaceData,Descriptor>::container);

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Domain Information Lattice = "+ box_string(lattice)+" Obstacle = "
					+box_string(obstacle_domain)+" and Wall = "+box_string(wall_domain);
				if(master){ std::cout << mesg << std::endl; }
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
				GuoOffLatticeModel3D<T,Descriptor>* model_ptr = model->clone();

				pcout << &model_ptr << std::endl;

				GetForceOnObjectFunctional3D<T,BoundaryType> force(model_ptr);

				//PlainReductiveBoxProcessingFunctional3D* func = nullptr;
				//func = dynamic_cast<PlainReductiveBoxProcessingFunctional3D*>(force.get());

				Box3D domain(0,0,0,0,0,0);
				domain = bc->getArg().getBoundingBox();
				pcout << "DOMAIN ["<<domain.x0 << "," <<domain.x1 <<"," <<domain.y0 << "," <<domain.y1<<","<<domain.z0<<","
					<<domain.z1<<"]"<<std::endl;

				std::vector<MultiBlock3D*> arg;
				arg.push_back(dynamic_cast<MultiBlock3D*>(model->generateOffLatticeInfo()));
				MultiBlock3D* bcarg = &bc->getArg();
				arg.push_back(bcarg);
				arg.push_back(dynamic_cast<MultiBlock3D*>(Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->clone()));
				/*
				std::map<plint,BlockLattice3D<T,Descriptor>*> blocks =
					Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBlockLattices();
				std::vector<AtomicBlock3D*> atomics;
				plint size = blocks.size();
				atomics.resize(size);
				atomics.reserve(size);
				atomics.push_back(dynamic_cast<AtomicBlock3D*>(model->generateOffLatticeInfo()));
				for(int i = 0; i<size; i++){
					atomics.push_back(dynamic_cast<AtomicBlock3D*>(blocks[i]));
				}
				pcout << atomics.size() << std::endl;
				pcout << &force << std::endl;
				pcout << &domain << std::endl;
				pcout << &atomics << std::endl;
				*/

				//force.processGenericBlocks(domain, atomics);

				applyProcessingFunctional(force, domain, arg);

				/*
				if(force != nullptr){ }
				else{ throw std::runtime_error("GetForceOnObjectFunctional3D not Properly Initialized"); }
				*/

				f += force.getForce();
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
				const T omega = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getOmega();
				const T rho_LB = (T)1.0;
				const T timeLB = Variables<T,BoundaryType,SurfaceData,Descriptor>::time;

				if(firstMove){ velocityFunc.initialize(mass, g, Constants<T>::obstacle.density, dt, dx, triangleSet); firstMove = false; }
				normalFunc.update(triangleSet);

				T factor = util::sqr(util::sqr(dx)) / util::sqr(dt);

				resetForceStatistics<T>(*Variables<T,BoundaryType,SurfaceData,Descriptor>::container);

				recomputeImmersedForce<T>(normalFunc, omega, rho_LB,
					*Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice,
					*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
					Constants<T>::envelopeWidth, Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBoundingBox(), true);

				Array<T,3> force = Array<T,3>(0,0,0);
				force = reduceImmersedForce<T>(*Variables<T,BoundaryType,SurfaceData,Descriptor>::container);

				T x = 0;
				T y = 0;
				T z = 0;
				for(int i = 0; i<numVertices; i++){
					Array<T,3> iVertex = triangleSet.getVertex(i);
					x += iVertex[0];
					y += iVertex[1];
					z += iVertex[2];
				}
				Array<T,3> center = Array<T,3>(x/numVertices, y/numVertices, z/numVertices);

				Array<T,3> torque = Array<T,3>(0,0,0);
				torque = reduceAxialTorqueImmersed(*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
										center, Array<T,3>(1,1,1));

				Array<T,3> ds = Array<T,3>(0,0,0);
				ds = velocityFunc.update(timeLB,force,torque,triangleSet);
				
				for (int i = 0; i < Constants<T>::ibIter; i++){
					indexedInamuroIteration<T>(velocityFunc,
									*Variables<T,BoundaryType,SurfaceData,Descriptor>::rhoBar,
									*Variables<T,BoundaryType,SurfaceData,Descriptor>::j,
									*Variables<T,BoundaryType,SurfaceData,Descriptor>::container,
									Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getTau(),
									true);
				}

				normalFunc.update(triangleSet);

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

