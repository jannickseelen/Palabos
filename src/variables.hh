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
				std::string mesg ="[DEBUG] Updating Paramters";
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
			else{ writeLogFile(p,"parameters");}
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
	std::unique_ptr<DEFscaledMesh<T> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createMesh(
		const ConnectedTriangleSet<T>& triangleSet, const plint& referenceDirection, const int& flowType)
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

			mesh.reset(new DEFscaledMesh<T>(*triangleSet.toTriangleSet(Constants<T>::precision), resolution, referenceDirection,
						Constants<T>::margin, Constants<T>::extraLayer));


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

			const MultiBlockManagement3D& management = voxelizedDomain.getMultiBlockManagement();
			BlockCommunicator3D* communicator = voxelizedDomain.getBlockCommunicator();
			CombinedStatistics* statistics = voxelizedDomain.getCombinedStatistics();
			MultiCellAccess3D<T,Descriptor>* access = defaultMultiBlockPolicy3D().getMultiCellAccess<T,Descriptor>();
			Dynamics<T,Descriptor>* d = dynamics.get();

			partial_lattice.reset(new MultiBlockLattice3D<T,Descriptor>(management, communicator, statistics, access, d));
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
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createBP()
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
		const VoxelizedDomain3D<T>& voxelizedDomain, const BoundaryProfiles3D<T,SurfaceData>& profile)
	{
		std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > flowShape(nullptr);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Triangle Flow Shape";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			flowShape.reset(
				new TriangleFlowShape3D<T,SurfaceData>(
					voxelizedDomain.getBoundary(),
					profile)
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
	std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > Variables<T,BoundaryType,SurfaceData,Descriptor>::createBC(
		GuoOffLatticeModel3D<T,Descriptor>* model, VoxelizedDomain3D<T>& voxelizedDomain, MultiBlockLattice3D<T,Descriptor>& lt)
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
					lt)
			);

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
	void Variables<T,BoundaryType,SurfaceData,Descriptor>::join()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Joining Lattices ";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("join").start();
			#endif

			Box3D fromDomain = Obstacle<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBoundingBox();
			Box3D toDomain = Wall<T,BoundaryType,SurfaceData,Descriptor>::lattice->getBoundingBox();

			lattice.reset(new MultiBlockLattice3D<T,Descriptor>(*Wall<T,BoundaryType,SurfaceData,Descriptor>::lattice));
			lattice->copyReceive(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::lattice,
				fromDomain, toDomain, modif::allVariables);
			defineDynamics(*lattice, lattice->getBoundingBox(), dynamics->clone());
			lattice->toggleInternalStatistics(false);

			rhoBar.reset(generateMultiScalarField<T>((MultiBlock3D&) *lattice, Constants<T>::envelopeWidth).release());
			rhoBar->toggleInternalStatistics(false);

			j.reset(generateMultiTensorField<T,3>((MultiBlock3D&) *lattice, Constants<T>::envelopeWidth).release());
			j->toggleInternalStatistics(false);

			rhoBarJarg.push_back(dynamic_cast<MultiBlock3D*>(lattice.get()));
			rhoBarJarg.push_back(dynamic_cast<MultiBlock3D*>(rhoBar.get()));
			rhoBarJarg.push_back(dynamic_cast<MultiBlock3D*>(j.get()));

			integrateProcessingFunctional(new ExternalRhoJcollideAndStream3D<T,Descriptor>(),lattice->getBoundingBox(), rhoBarJarg, 0);
			integrateProcessingFunctional(new BoxRhoBarJfunctional3D<T,Descriptor>(), lattice->getBoundingBox(), rhoBarJarg, 2);

			lattice->periodicity().toggleAll(false);
			rhoBar->periodicity().toggleAll(false);
			j->periodicity().toggleAll(false);

			container.reset(new MultiContainerBlock3D(*rhoBar));

			Wall<T,BoundaryType,SurfaceData,Descriptor>::bc->insert(rhoBarJarg);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bc->insert(rhoBarJarg);

			initializeAtEquilibrium(*lattice, lattice->getBoundingBox(), (T)1.0, Array<T,3>((T) 0, (T) 0, (T) 0));

			T latticeU = p.getLatticeU();
			T Re = p.getRe();
			plint resolution = p.getResolution();
			Box3D domain = lattice->getBoundingBox();
			T lx = lattice->getNx();
			T ly = lattice->getNy();
			T lz = lattice->getNz();
			p = IncomprFlowParam<T>(latticeU, Re, resolution, lx, ly, lz);
			writeLogFile(p, "parameters");

			#ifdef PLB_DEBUG
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

			Wall<T,BoundaryType,SurfaceData,Descriptor>::mesh = createMesh(Wall<T,BoundaryType,SurfaceData,Descriptor>::triangleSet,
				Constants<T>::wall.referenceDirection, Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::tb = createTB(*Wall<T,BoundaryType,SurfaceData,Descriptor>::mesh);
			Wall<T,BoundaryType,SurfaceData,Descriptor>::location = Wall<T,BoundaryType,SurfaceData,Descriptor>::tb->getPhysicalLocation();

			Wall<T,BoundaryType,SurfaceData,Descriptor>::vd = createVoxels(*Wall<T,BoundaryType,SurfaceData,Descriptor>::tb,
				Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType, Constants<T>::wall.dynamicMesh);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::lattice = createLattice(*Wall<T,BoundaryType,SurfaceData,Descriptor>::vd);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::bp = createBP();

			Wall<T,BoundaryType,SurfaceData,Descriptor>::fs = createFS(*Wall<T,BoundaryType,SurfaceData,Descriptor>::vd,
				*Wall<T,BoundaryType,SurfaceData,Descriptor>::bp);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::model = createModel(Wall<T,BoundaryType,SurfaceData,Descriptor>::fs.get(),
				Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Wall<T,BoundaryType,SurfaceData,Descriptor>::bc = createBC(Wall<T,BoundaryType,SurfaceData,Descriptor>::model.get(),
				*Wall<T,BoundaryType,SurfaceData,Descriptor>::vd, *Wall<T,BoundaryType,SurfaceData,Descriptor>::lattice);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::mesh = createMesh(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::triangleSet,
				Constants<T>::obstacle.referenceDirection, Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb = createTB(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::mesh);
			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::location = Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb->getPhysicalLocation();

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd = createVoxels(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb,
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType, Constants<T>::obstacle.dynamicMesh);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::lattice = createLattice(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bp = createBP();

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::fs = createFS(*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd,
				*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bp);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::model = createModel(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::fs.get(),
				Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType);

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bc = createBC(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::model.get(),
				*Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd, *Obstacle<T,BoundaryType,SurfaceData,Descriptor>::lattice);

			join();

			Obstacle<T,BoundaryType,SurfaceData,Descriptor>::moveToStart();

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

			MultiScalarField3D<T> r = *computeDensity(*lattice);
			density.push_back(r);

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

			save();

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


