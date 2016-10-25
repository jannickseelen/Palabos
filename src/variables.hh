#ifndef VARIABLES_HH
#define VARIABLES_HH

#include "variables.h"
#include <palabos3D.hh>
#include "myheaders3D.hh"

#include <thread>
#include <chrono>

namespace plb{

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Variables<T,BoundaryType,SurfaceData,Descriptor>::Variables()
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing Variables";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->v.reset(this);
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Variables already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Variables<T,BoundaryType,SurfaceData,Descriptor>::~Variables()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Variables";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		objCount--;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::initialize()
	{
		try{
			resolution = 0; gridLevel=0; reynolds=0;
			location = Array<T,3>();
			dx = 1;
			dt = 1;
			iter = 0;
			nprocs = 0;
			nprocs_side = 0;
			master = global::mpi().isMainProcessor();
			nprocs = global::mpi().getSize();
			nprocs_side = (int)cbrt(nprocs);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::update(const plint& _gridLevel, const plint& _reynolds)
	{
		try{
			Constants<T>* c = Constants<T>::c.get();
			#ifdef PLB_DEBUG
				std::string mesg ="[DEBUG] Updating Parameters";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			gridLevel = _gridLevel;
			resolution = Constants<T>::physical.resolution * util::twoToThePowerPlint(_gridLevel);
			scaled_u0lb = Constants<T>::lb.u * util::twoToThePowerPlint(_gridLevel);
			reynolds = _reynolds;
			p = IncomprFlowParam<T>(Constants<T>::physical.u,scaled_u0lb,reynolds,Constants<T>::physical.length,
				Constants<T>::physical.resolution,Constants<T>::lb.lx,Constants<T>::lb.ly,Constants<T>::lb.lz);
			dynamics.reset(new IncBGKdynamics<T,Descriptor>(p.getOmega()));
			dx = p.getDeltaX();
			dt = p.getDeltaT();
			T U = p.getLatticeU();
			scalingFactor = (T)(resolution)/dx;
			//for(int i = 0; i<rhoBarJarg.size(); i++){ delete rhoBarJarg[i];}
			rhoBarJarg.clear();
			velocity.clear();
			vorticity.clear();
			density.clear();
			//lattice.reset(nullptr);
			//rhoBar.reset(nullptr);
			//j.reset(nullptr);
			//container.reset(nullptr);
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Reynolds="+std::to_string(reynolds);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Grid Level="+std::to_string(gridLevel);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Resolution="+std::to_string(resolution);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Lattice U="+std::to_string(U);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Lattice dt="+std::to_string(dt);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Lattice dx="+std::to_string(dx);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Scaling Factor="+std::to_string(scalingFactor);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Done Updating Parameters";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			if(U == 0 || dx == 0 || dt == 0){
				writeLogFile(p, "WRONG PARAMATERS");
				throw std::runtime_error("InComprFlowParam not set Correctly");
			}
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	bool Variables<T,BoundaryType,SurfaceData,Descriptor>::checkConvergence()
	{
		try{
			util::ValueTracer<T> tracer(p.getLatticeU(), p.getDeltaX(), Constants<T>::epsilon);
			return tracer.hasConverged();
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return false;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	T Variables<T,BoundaryType,SurfaceData,Descriptor>::getRho(const T& temp)
	{
		T rho_lb = (T)1.0;
		try{
			const T dx = p.getDeltaX();
			T rho = 2579.3 - 0.6237*temp;
			rho_lb = rho*dx*dx*dx;
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return rho_lb;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<DEFscaledMesh<T> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createMesh(
		ConnectedTriangleSet<T>& triangleSet, const plint& referenceDirection, const int& flowType)
	{
		std::unique_ptr<DEFscaledMesh<T> > mesh(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Mesh";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").start();
			#endif
			// Create Mesh

			T dx = p.getDeltaX();
			T alpha = 1/dx;
			TriangleSet<T> normalTriangleSet = *triangleSet.toTriangleSet(Constants<T>::precision);

			Cuboid<T> cube = normalTriangleSet.getBoundingCuboid();
			Array<T,3> lowerLeftCorner = cube.lowerLeftCorner;
			Array<T,3> upperRightCorner = cube.upperRightCorner;

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Bounded Cuboid BEFORE Scaling Lower Left Corner "+array_string(lowerLeftCorner)+" Upper Right Corner "+
					array_string(upperRightCorner)+" in physical units";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			normalTriangleSet.scale(alpha);

			cube = normalTriangleSet.getBoundingCuboid();
			lowerLeftCorner = cube.lowerLeftCorner;
			upperRightCorner = cube.upperRightCorner;
			nx = upperRightCorner[0]+1;
			ny = upperRightCorner[1]+1;
			nz = upperRightCorner[2]+1;

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Bounded Cuboid AFTER Scaling Lower Left Corner "+array_string(lowerLeftCorner)+" Upper Right Corner "+
					array_string(upperRightCorner)+" in dimensionless units";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			mesh.reset(new DEFscaledMesh<T>(normalTriangleSet, resolution, referenceDirection,
						Constants<T>::margin, Constants<T>::extraLayer));

			triangleSet = ConnectedTriangleSet<T>(normalTriangleSet);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Mesh address= "+adr_string(mesh.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return mesh;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<TriangleBoundary3D<T> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createTB(const DEFscaledMesh<T>& mesh)
	{
		std::unique_ptr<TriangleBoundary3D<T> > triangleBoundary(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Triangle Boundary";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
			// Create Mesh
			triangleBoundary.reset(new TriangleBoundary3D<T>(mesh,true));
			triangleBoundary->getMesh().inflate();
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] TriangleBoundary address= "+adr_string(triangleBoundary.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return triangleBoundary;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<VoxelizedDomain3D<T> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createVoxels(const TriangleBoundary3D<T>& tb,
		const int& flowType, const bool& dynamic)
	{
		std::unique_ptr<VoxelizedDomain3D<T> > voxelizedDomain(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Voxelized Domain";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				std::cout << "[DEBUG] TB Address=" << &tb << " FlowType="<<flowType<<" ExtraLayer="<<Constants<T>::extraLayer <<
					" Borderwidth= "<<Constants<T>::borderWidth<<" EnvelopeWidth= "<<Constants<T>::envelopeWidth<<" Blocksize= "<<
					Constants<T>::blockSize<<" GridLevel= "<<gridLevel<<" Dynamic Mesh= "<<dynamic<<std::endl;
				global::timer("boundary").restart();
			#endif
			voxelizedDomain.reset(
				new VoxelizedDomain3D<T>(
					tb,
					flowType,
					Constants<T>::extraLayer,
					Constants<T>::borderWidth,
					Constants<T>::envelopeWidth,
					Constants<T>::blockSize,
					gridLevel,
					dynamic));
			MultiScalarField3D<int> flagMatrix((MultiBlock3D&)voxelizedDomain->getVoxelMatrix());
			if(flowType == voxelFlag::inside){
				setToConstant(flagMatrix, voxelizedDomain->getVoxelMatrix(),	voxelFlag::inside, flagMatrix.getBoundingBox(), 1);
				setToConstant(flagMatrix, voxelizedDomain->getVoxelMatrix(), voxelFlag::innerBorder, flagMatrix.getBoundingBox(), 1);
			}
			else{
				setToConstant(flagMatrix, voxelizedDomain->getVoxelMatrix(),	voxelFlag::outside, flagMatrix.getBoundingBox(), 1);
				setToConstant(flagMatrix, voxelizedDomain->getVoxelMatrix(), voxelFlag::outerBorder, flagMatrix.getBoundingBox(), 1);
			}
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] VoxelizedDomain  address= "+adr_string(voxelizedDomain.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return voxelizedDomain;
	}

	/*
	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createLattice(
		VoxelizedDomain3D<T>& voxelizedDomain)
	{
		std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > partial_lattice(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Partial Lattice";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			Dynamics<T,Descriptor>* d = dynamics.get();

			partial_lattice.reset(new MultiBlockLattice3D<T,Descriptor>(nx, ny, nz, d));
			partial_lattice->toggleInternalStatistics(false);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Partial Lattice address= "+adr_string(partial_lattice.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return partial_lattice;
	}*/

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createBP(
		const TriangleBoundary3D<T>& tb)
	{
		std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > profile(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Boundary Profile";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			profile.reset(new BoundaryProfiles3D<T,SurfaceData>());
			profile->setWallProfile(NoSlipProfile3D<T>().clone());
			plint numTriangles = tb.getMesh().getNumTriangles();
			int n=(int)numTriangles;
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Defining Profiles for "+safe_string(n)+" surface triangles";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			for(plint t = 0; t<numTriangles; t++){
				plint tag = tb.getTag(t);
				profile->defineProfile(tag, NoSlipProfile3D<T>().clone());
			}

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Boundary Profile address= "+adr_string(profile.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return profile;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createFS(
		const VoxelizedDomain3D<T>& voxelizedDomain, const BoundaryProfiles3D<T,SurfaceData>* profile)
	{
		std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > flowShape(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Triangle Flow Shape";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			if(profile == nullptr){ throw std::runtime_error("Boundary Profiles where destroyed!"); }

			flowShape.reset(
				new TriangleFlowShape3D<T,SurfaceData>(
					voxelizedDomain.getBoundary(),
					*profile)
			);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Triangle Flow Shape address= "+adr_string(flowShape.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return flowShape;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<GuoOffLatticeModel3D<T,Descriptor> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createModel(
		TriangleFlowShape3D<T,SurfaceData>* flowShape, const int& flowType)
	{
		std::unique_ptr<GuoOffLatticeModel3D<T,Descriptor> > model(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Model";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			model.reset(
				new GuoOffLatticeModel3D<T,Descriptor>(
					flowShape,
					flowType)
			);

			model->setVelIsJ(true); // When the incompressible BGK model is used, velocity equals momentum.
			model->selectUseRegularizedModel(true);
			model->selectComputeStat(false);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Model address= "+adr_string(model.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return model;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::createLattice(const VoxelizedDomain3D<T>& wallVoxels,
		const VoxelizedDomain3D<T>& obstacleVoxels)
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Joining Lattices ";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("join").start();
			#endif

			MultiScalarField3D<int> wallMatrix = wallVoxels.getVoxelMatrix();
			MultiScalarField3D<int> obstacleMatrix = obstacleVoxels.getVoxelMatrix();

			Box3D fromDomain = Obstacle<T,BoundaryType,SurfaceData,Descriptor>::getDomain();
			Box3D toDomain = Wall<T,BoundaryType,SurfaceData,Descriptor>::getDomain();

			wallMatrix.copyReceive(obstacleMatrix,fromDomain,toDomain,modif::allVariables);

			lattice.reset(generateMultiBlockLattice<T,Descriptor>(wallMatrix, Constants<T>::envelopeWidth, dynamics->clone()).release());

			T resolution = Constants<T>::physical.resolution * util::twoToThePowerPlint(gridLevel);
			T scaled_u0lb = Constants<T>::lb.u / util::twoToThePowerPlint(gridLevel);
			p = IncomprFlowParam<T>(Constants<T>::physical.u, scaled_u0lb, reynolds, Constants<T>::physical.length,
									resolution, lattice->getNx(), lattice->getNy(), lattice->getNz());
			std::string fileName = "parameters_Re="+std::to_string((int)reynolds)+"_GridLvL="+std::to_string((int)gridLevel);
			writeLogFile(p, fileName);

			defineDynamics(*lattice, lattice->getBoundingBox(), dynamics->clone());
			lattice->toggleInternalStatistics(false);
			defineDynamics(*lattice, wallMatrix, lattice->getBoundingBox(), new NoDynamics<T,Descriptor>, voxelFlag::outside);
			defineDynamics(*lattice, obstacleMatrix, lattice->getBoundingBox(), new NoDynamics<T,Descriptor>, voxelFlag::inside);

			rhoBar.reset(generateMultiScalarField<T>((MultiBlock3D&) *lattice, Constants<T>::envelopeWidth).release());
			rhoBar->toggleInternalStatistics(false);

			j.reset(generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, Constants<T>::envelopeWidth).release());
			j->toggleInternalStatistics(false);

			rhoBarJarg.clear();
			rhoBarJarg.push_back(dynamic_cast<MultiBlock3D*>(lattice.get()));
			rhoBarJarg.push_back(dynamic_cast<MultiBlock3D*>(rhoBar.get()));
			rhoBarJarg.push_back(dynamic_cast<MultiBlock3D*>(j.get()));

			integrateProcessingFunctional(new ExternalRhoJcollideAndStream3D<T,Descriptor>(),lattice->getBoundingBox(), rhoBarJarg, 0);
			integrateProcessingFunctional(new BoxRhoBarJfunctional3D<T,Descriptor>(), lattice->getBoundingBox(), rhoBarJarg, 3);

			lattice->periodicity().toggleAll(false);
			rhoBar->periodicity().toggleAll(false);
			j->periodicity().toggleAll(false);

			container = new MultiContainerBlock3D(*rhoBar);

			T lx = lattice->getNx();
			T ly = lattice->getNy();
			T lz = lattice->getNz();

			Box3D domain = lattice->getBoundingBox();

			TriangularSurfaceMesh<T> mesh = Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb->getMesh();
			//ConnectedTriangleSet<T> triangles = ConnectedTriangleSet<T>(mesh.toTriangleSet(Constants<T>::precision));

			T numVertices = Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb->getMesh().getNumVertices();
			//T numVertices = triangles.getNumVertices();
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::numVertices = numVertices;

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vertices.clear();
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vertices.resize(numVertices);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vertices.reserve(numVertices);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::unitNormals.clear();
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::unitNormals.resize(numVertices);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::unitNormals.reserve(numVertices);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::areas.clear();
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::areas.resize(numVertices);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::areas.reserve(numVertices);

			const bool weightedArea = false;

			for(int i = 0; i < numVertices; i++){
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vertices[i] = mesh.getVertex(i);
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::areas[i] = mesh.computeVertexArea(i);
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::unitNormals[i] = mesh.computeVertexNormal(i,weightedArea);
				/*
				Array<T,3> v = Array<T,3>(0,0,0);
				v = triangles.getVertex(i);
				Array<T,3> n = Array<T,3>(0,0,0);
				T a = 0;
				triangles.computeVertexAreaAndUnitNormal(i, a, n);
				vertices[i] = v;
				areas[i] = a;
				unitNormals[i] = n;
				*/
			}

			// Integrate the immersed boundary processors in the lattice multi-block.
			std::vector<MultiBlock3D*> args;
			plint pl = 4;

			args.resize(0);
			args.push_back(container);
			integrateProcessingFunctional(new InstantiateImmersedWallData3D<T>(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vertices,
											Obstacle<T,BoundaryType,SurfaceData,Descriptor>::areas,
											Obstacle<T,BoundaryType,SurfaceData,Descriptor>::unitNormals),
											container->getBoundingBox(), *lattice, args, pl);
			//lattice->executeInternalProcessors(pl);
			//instantiateImmersedWallData(vertices, areas, *container);

			pl++;
			// Update the Velocity Function once
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::velocityFunc.update(p,(T)0,Array<T,3>(0,0,0),Array<T,3>(0,0,0),
					Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb.get(), lattice->getBoundingBox());

			for (plint i = 0; i < Constants<T>::ibIter; i++) {
				args.resize(0);
				args.push_back(rhoBar.get());
				args.push_back(j.get());
				args.push_back(container);
				integrateProcessingFunctional(
					new IndexedInamuroIteration3D<T,SurfaceVelocity<T> >(
						Obstacle<T,BoundaryType,SurfaceData,Descriptor>::velocityFunc, p.getTau(), true),
						rhoBar->getBoundingBox(), *lattice, args, pl);
				pl++;
			}
			Box3D newDomain = lattice->getBoundingBox();

			if(newDomain.x0 > fromDomain.x0 || newDomain.x1 < fromDomain.x1
			|| newDomain.y0 > fromDomain.y0 || newDomain.y1 < fromDomain.y1
			|| newDomain.z0 > fromDomain.z0 || newDomain.z1 < fromDomain.z1
			|| newDomain.x0 > toDomain.x0 || newDomain.x1 < toDomain.x1
			|| newDomain.y0 > toDomain.y0 || newDomain.y1 < toDomain.y1
			|| newDomain.z0 > toDomain.z0 || newDomain.z1 < toDomain.z1)
			{
				std::string ex = "[ERROR] Domain mismatch Lattice = "+ box_string(newDomain)+" but obstacle = "
				+box_string(fromDomain)+" and wall = "+box_string(toDomain);
				throw std::runtime_error(ex);
			}

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Domain Information Lattice = "+ box_string(newDomain)+" Obstacle = "
					+box_string(fromDomain)+" and Wall = "+box_string(toDomain);
				if(master){ std::cout << mesg << std::endl; }
				global::log(mesg);
				mesg = "[DEBUG] Domain= "+ box_string(domain)+" Nx= "+std::to_string(lx)+" Ny= "+std::to_string(ly)+" Nz= "+std::to_string(lz);
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Done Joining Lattices time="+std::to_string(global::timer("join").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("join").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createBC(
		GuoOffLatticeModel3D<T,Descriptor>* model, VoxelizedDomain3D<T>& voxelizedDomain)
	{
		std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > boundaryCondition(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating BoundaryCondition";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			boundaryCondition.reset(
					new OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>(
					model,
					voxelizedDomain,
					*lattice)
			);

			boundaryCondition->insert(rhoBarJarg);

			#ifdef PLB_DEBUG
				mesg ="[DEBUG] BoundaryCondition address= "+adr_string(boundaryCondition.get());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return boundaryCondition;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::initializeLattice()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Initializing Lattice";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("ini").start();
			#endif

			T iniT = Constants<T>::initialTemperature;
			T rho_lb = getRho(iniT);
			Array<T,3> iniV = Array<T,3>(0,0,0);

			initializeAtThermalEquilibrium(*lattice, lattice->getBoundingBox(), rho_lb, iniV, iniT);
			//lattice->initialize();
			applyProcessingFunctional(new BoxRhoBarJfunctional3D<T,Descriptor>(),lattice->getBoundingBox(), rhoBarJarg);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::bc->apply(rhoBarJarg);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bc->apply(rhoBarJarg);

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Initializing Lattice time="+std::to_string(global::timer("join").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("ini").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::makeParallel()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Parallelizing Lattice ";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("parallel").start();
			#endif
			Box3D box = lattice->getBoundingBox();
			std::map<plint, BlockLattice3D< T, Descriptor > * > blockMap = lattice->getBlockLattices();
			plint size = blockMap.size();
			std::vector< std::vector<Box3D> > domains;
			domains.resize(size);

			for(int i=0; i<size; i++){
				std::vector<Box3D> block;
				plint nx = blockMap[i]->getNx()-1;
				plint ny = blockMap[i]->getNy()-1;
				plint nz = blockMap[i]->getNz()-1;
				for(int x = 0; x<nx; x++){
					for(int y = 0; y<nx; y++){
						for(int z=0; z<nz; z++){
							Box3D cell(x,x+1,y,y+1,z,z+1);
							block.push_back(cell);
						}
					}
				}
				domains.push_back(block);
			}
			ParallellizeByCubes3D parallel(domains, box, nprocs_side, nprocs_side, nprocs_side);
			parallel.parallelize();
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Parallelizing Lattice time="+std::to_string(global::timer("parallel").getTime());
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("parallel").stop();
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::setLattice()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing Main Lattice";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			std::unique_ptr<DEFscaledMesh<T> > wall_mesh = createMesh(Wall<T,BoundaryType,SurfaceData,Descriptor>::triangleSet,
				Constants<T>::wall.referenceDirection, Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			std::unique_ptr<DEFscaledMesh<T> > obstacle_mesh = createMesh(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::triangleSet,
				Constants<T>::obstacle.referenceDirection, Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::tb = createTB(*wall_mesh);
			Wall<T,BoundaryType,SurfaceData,Descriptor>::location = Wall<T,BoundaryType,SurfaceData,Descriptor>::tb->getPhysicalLocation();

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb = createTB(*obstacle_mesh);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::location = Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb->getPhysicalLocation();

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::moveToStart();

			Wall<T,BoundaryType,SurfaceData,Descriptor>::vd = createVoxels(*Wall<T,BoundaryType,SurfaceData,Descriptor>::tb,
				Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType, Constants<T>::wall.dynamicMesh);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd = createVoxels(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb,
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType, Constants<T>::obstacle.dynamicMesh);

			createLattice(*Wall<T,BoundaryType,SurfaceData,Descriptor>::vd, *Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::bp = createBP(*Wall<T,BoundaryType,SurfaceData,Descriptor>::tb);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::fs = createFS(*Wall<T,BoundaryType,SurfaceData,Descriptor>::vd,
				Wall<T,BoundaryType,SurfaceData,Descriptor>::bp.get());

			Wall<T,BoundaryType,SurfaceData,Descriptor>::model = createModel(Wall<T,BoundaryType,SurfaceData,Descriptor>::fs.get(),
				Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::bc = createBC(Wall<T,BoundaryType,SurfaceData,Descriptor>::model.get(),
				*Wall<T,BoundaryType,SurfaceData,Descriptor>::vd);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bp = createBP(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::fs = createFS(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd,
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bp.get());

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::model = createModel(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::fs.get(),
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bc = createBC(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::model.get(),
				*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd);

			initializeLattice();

			//makeParallel();

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Constructing Main Lattice";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::save()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Saving Data";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			lattice->toggleInternalStatistics(false);

			//MultiTensorField3D<T, 3> v = *computeVelocity(*lattice);
			//velocity.push_back(v);

			//MultiTensorField3D<T, 3> w = *computeVorticity(v);
			//vorticity.push_back(w);

			//MultiScalarField3D<T> r = *computeDensity(*lattice);
			//density.push_back(r);

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Saving Data";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::updateLattice()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating Main Lattice";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			lattice->toggleInternalStatistics(false);

			//save();

			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Updating Main Lattice";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}



} // namespace plb

#endif // VARIABLES_HH


