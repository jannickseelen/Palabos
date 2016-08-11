//
//  test.cpp
//  Palabos
//
//  Created by MAC on 10/02/16.
//  Copyright © 2016 TUDelft. All rights reserved.
//
//
// C++ INCLUDES
#include <execinfo.h> // For complete stacktrace on exception see GNU libc manual
#include <iostream>
#include <math.h>
#include <memory>
#include <thread>
#include <chrono>
#include <exception>
#include <string>
#include <sstream>
// PALABOS INCLUDES
#include <offLattice/offLatticeBoundaryCondition3D.hh>
#include <offLattice/triangleSet.hh>
#include <offLattice/guoOffLatticeModel3D.hh>
#include <basicDynamics/isoThermalDynamics.hh>
#include <latticeBoltzmann/nearestNeighborLattices3D.hh>
#include <dataProcessors/dataAnalysisFunctional3D.hh>
#include <dataProcessors/dataInitializerFunctional3D.hh>
#include <multiBlock/multiBlockLattice3D.hh>
#include <multiBlock/multiDataProcessorWrapper3D.hh>
#include <multiGrid/parallelizer3D.h>
#include <atomicBlock/dataProcessingFunctional3D.hh>
#include <core/blockLatticeBase3D.hh>
#include <core/units.h>
#include <core/plbInit.hh>
#include <core/dynamicsIdentifiers.hh>
#include <core/plbLogFiles.h>
#include <libraryInterfaces/TINYXML_xmlIO.hh>
#include <io/vtkDataOutput.hh>
#include <io/imageWriter.hh>
#include <algorithm/benchmarkUtil.hh>
#include <modules/mpiDataManager.hh>
#include <modules/LBMexceptions.hh>
#include <modules/debug.hh>

#define Descriptor plb::descriptors::D3Q27Descriptor

#ifdef PLB_DEBUG
	#define NDEBUG
#endif

namespace plb{

template<typename T>
class Point{
// Default Constructor
public:
	Point(){
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}
	Point(const T &_x, const T &_y, const T &_z):x(_x),y(_y),z(_z){}
	Point(const Array<T,3> &array):x(array[0]),y(array[1]),z(array[2]){}
	Point(const Point &p):x(p.x),y(p.y),z(p.z){}
// Default Destructor
	~Point(){}
// Methods
	T distance(const Point &p){
		T d = 0;
		d = sqrt(pow((this->x - p.x),2) + pow((this->y - p.y),2) + pow((this->z - p.z),2));
		return d;
	}
// Properties
T x;
T y;
T z;
};

template<typename T>
class Triangle{
// Default Constructor
public:
	Triangle(){} // Default constructor that calls the default point constructor for a,b and c
	Triangle(const Point<T> &_a, const Point<T> &_b, const Point<T> &_c):a(_a),b(_b),c(_c){}
	Triangle(const Array<Array<T,3>,3> &array){
			this->a = Point<T>(array[0]);
			this->b = Point<T>(array[1]);
			this->c = Point<T>(array[2]);
	}
// Default Destructor
	~Triangle(){};
// Methods
	Point<T> getCentroid(){
		Point<T> center(1/3*(this->a.x + this->b.x + this->c.x), 1/3*(this->a.y + this->b.y + this->c.y), 1/3*(this->a.z + this->b.z + this->c.z));
		return c;
	}
	double area() { // Heron's formula
		double area = 0; double AB=0; double BC=0; double CA=0; double s = 0;
		AB = this->a.distance(this->b);
		BC = this->b.distance(this->c);
		CA = this->c.distance(this->a);
		s = (AB + BC + CA)/2;
		area = sqrt(s*(s-AB)*(s-BC)*(s-CA));
		return area;
	}
// Properties
	Point<T> a, b, c;
};

template<typename T>
class Pyramid{
// Default Constructor
public:
	Pyramid(){} // Default constructor that calls the default point constructor for a,b and c
	Pyramid(const Triangle<T> &b, const Point<T> &t):base(b),top(t){}
// Default Destructor
	~Pyramid(){};
// Methods
	T volume() {
		T volume = 0; T area = 0; T height = 0;
		Point<T> centroid = base.getCentroid();
		area = base.area();
		height = centroid.distance(this->top);
		volume = (1/3)*area*height;
		return volume;
	}
// Properties
	Triangle<T> base;
	Point<T> top;
};

template<typename T, class BoundaryType>
class Constants{
private:
	static Constants<T,BoundaryType>* c;
	Constants(){ //Private Constructor for Static Instance
		this->parameterXmlFileName = ""; this->u0lb=0; this->epsilon=0; this->minRe = 0; this->maxRe=0;this->maxGridLevel=0;
		this->referenceResolution=0;
		this->margin=0; this->borderWidth=0; this->extraLayer=0; this->blockSize=0; this->envelopeWidth=0; this->initialTemperature=0;
		this->gravitationalAcceleration=0; this->dynamicMesh = false; this->master = false;
		this->test = true; this->testRe = 0; this->testTime = 20; this->maxT = 0; this->imageSave = 0; this->testIter = 0;
		this->master = global::mpi().isMainProcessor();
	}
	// Default Destructor
	~Constants(){
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Constants DESTRUCTOR was called";
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
			throw std::runtime_error("Constants Destructor was Called");
		#endif
		delete c;
	}
public:
// Methods
	static Constants* getInstance(){ if(!c || c == nullptr){ c = new Constants<T,BoundaryType>();} return c; }

	void initialize(const std::string& fileName){
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Constants";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->parameterXmlFileName = fileName;
			XMLreader r(fileName);
			std::string log_dir;
			r["dir"]["log"].read(log_dir);
			global::directories().setLogOutDir(log_dir);
			r["lbm"]["u0lb"].read(this->u0lb);
			double x = 1;
			double y = 3;
			double cs = sqrt(x / y);
			double mach = this->u0lb / cs;
			double maxMach = 0.1;
			if(mach > maxMach){throw machEx;}
			r["lbm"]["maxRe"].read(this->maxRe);
			r["lbm"]["minRe"].read(this->minRe);
			r["lbm"]["epsilon"].read(this->epsilon);
			r["simulation"]["initialTemperature"].read(this->initialTemperature);
			r["refinement"]["margin"].read(this->margin);
			r["refinement"]["borderWidth"].read(this->borderWidth);
			if(this->margin < this->borderWidth){throw marginEx;} // Requirement: margin>=borderWidth.
			r["refinement"]["maxGridLevel"].read(this->maxGridLevel);
			r["refinement"]["extraLayer"].read(this->extraLayer);
			r["refinement"]["blockSize"].read(this->blockSize);
			r["refinement"]["envelopeWidth"].read(this->envelopeWidth);
			r["refinement"]["dynamicMesh"].read(this->dynamicMesh);
			r["refinement"]["referenceResolution"].read(this->referenceResolution);
			int tmp = 0;
			r["simulation"]["test"].read(tmp);
			if(tmp == 1){ this->test = true; }else{ this->test == false;}
			r["simulation"]["testRe"].read(this->testRe);
			r["simulation"]["testTime"].read(this->testTime);
			r["simulation"]["maxT"].read(this->maxT);
			r["simulation"]["imageSave"].read(this->imageSave);
			r["simulation"]["testIter"].read(this->testIter);
			if(this->test){ this->minRe = this->testRe; this->maxRe = this->testRe+1;  }
			double ratio = 3;
			// Fill the 2D array with standard values
			for(plint grid = 0; grid <= this->maxGridLevel; grid++){
				T resolution = this->referenceResolution * util::twoToThePowerPlint(grid);
				T scaled_u0lb = this->u0lb * util::twoToThePowerPlint(grid);
				mach = scaled_u0lb / cs;
				if(mach > maxMach){std::cout<<"Local Mach= "<<mach<<"\n"; throw localMachEx;}
				if(resolution == 0){throw resolEx;}
				if(this->test){
					IncomprFlowParam<T> p = IncomprFlowParam<T>(scaled_u0lb,testRe,resolution,1,1,1);
					// Check local speed of sound constraint
					T dt = p.getDeltaT();
					T dx = p.getDeltaX();
					if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
				}
				else{
					for(int reynolds = 0; reynolds <= this->maxRe; reynolds++){
						IncomprFlowParam<T> p = IncomprFlowParam<T>(scaled_u0lb,reynolds,resolution,1,1,1);
						// Check local speed of sound constraint
						T dt = p.getDeltaT();
						T dx = p.getDeltaX();
						if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
					}
				}
			}
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Creating Constants";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}
// Methods
// Properties
	std::string parameterXmlFileName;
	T u0lb, epsilon, maxT, imageSave;
	plint testIter, testRe, testTime, maxRe, minRe, maxGridLevel, referenceResolution, margin;
	plint borderWidth, extraLayer, blockSize, envelopeWidth;
	T initialTemperature, gravitationalAcceleration;
	bool dynamicMesh, test;
private:
	bool master;
};
template<typename T, class BoundaryType>
Constants<T,BoundaryType> *Constants<T,BoundaryType>::c=0;

template<typename T, class BoundaryType>
class Wall{
private:
	static Wall<T, BoundaryType>* w;
	Wall(){// Private Constructor for static Instance
		this->referenceDirection = 0;
		this->flowType = 0;
		this->temperature = 0; this->density=0;
		this->boundaryCondition = nullptr;
		this->c = Constants<T,BoundaryType>::getInstance();
		this->master = global::mpi().isMainProcessor();
	}
	~Wall(){
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Wall DESTRUCTOR was called";
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		throw std::runtime_error("Wall Destructor was Called");
		delete c;
		delete boundaryCondition;
	}
public:
// Methods
	static Wall* getInstance(){ if(!w || w == nullptr){ w = new Wall<T,BoundaryType>();} return w; }

	void initialize(){
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Initializing Wall";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			std::string meshFileName, material;
			int i = 0;
			XMLreader r(this->c->parameterXmlFileName);
			r["wall"]["meshFileName"].read(meshFileName);
			r["wall"]["referenceDirection"].read(this->referenceDirection);
			r["wall"]["initialTemperature"].read(this->temperature);
			r["wall"]["material"].read(material);
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] MeshFileName ="+meshFileName;
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			#ifdef PLB_MPI_PARALLEL
				if(global::mpi().isMainProcessor()){
					this->triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
					global::mpiData().sendTriangleSet<T>(this->triangleSet);}
				else{ this->triangleSet = global::mpiData().receiveTriangleSet<T>(); }
			#else
				this->triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
			#endif
			this->flowType = voxelFlag::inside;
			this->temperature = this->c->initialTemperature;
			if(material.compare("AL203")==0){ i = 1; }
			switch(i){
				case(0):	throw "Wall Material not Properly Defined! Modify input parameters.xml";
				case(1):	this->density = 3840; //[kg/m3];
			}
			#ifdef PLB_DEBUG
				mesg="[DEBUG] Number of triangles in Mesh = "+std::to_string(this->triangleSet.getTriangles().size());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg="[DEBUG] Done Initializing Wall";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}
	void setMesh(DEFscaledMesh<T>* fromMesh){
		this->mesh = fromMesh;
	}

	void setBoundary(OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* fromBoundary){
		this->boundaryCondition = fromBoundary;
	}
// Attributes
	plint referenceDirection;
	int flowType;
	T temperature, density;
	TriangleSet<T> triangleSet;
	DEFscaledMesh<T>* mesh;
	OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* boundaryCondition;
private:
	bool master;
	const Constants<T,BoundaryType>* c;
};
template<typename T, class BoundaryType>
Wall<T, BoundaryType> *Wall<T, BoundaryType>::w=0;

template<typename T, class BoundaryType>
class Obstacle{
private:
	static Obstacle<T, BoundaryType>* o;
	Obstacle(){ // Constructor for static Instance
		this->referenceDirection = 0;
		this->flowType = 0;
		this->density=0; this->mass=0; this->volume=0; this->temperature=0;
		this->center = Point<T>(); this->position = Point<T>();
		this->rotation = Array<T,3>();
		this->velocity = Array<T,3>();
		this->rotationalVelocity = Array<T,3>();
		this->acceleration = Array<T,3>();
		this->rotationalAcceleration = Array<T,3>();
		this->boundaryCondition = nullptr;
		this->c = Constants<T,BoundaryType>::getInstance();
		this->master = global::mpi().isMainProcessor();
	}
	// Default Destructor
	~Obstacle(){
		#ifdef PLB_DEBUG
			std::string mesg ="[DEBUG] Obstacle DESTRUCTOR was called";
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
			throw std::runtime_error("Obstacle Destructor was Called");
		#endif
		delete c;
		delete boundaryCondition;
		delete o;
	}
public:
	static Obstacle* getInstance(){ if(!o || o == nullptr){ o = new Obstacle<T,BoundaryType>(); } return o; }

	void initialize(){
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Initializing Obstacle";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			T x=0; T y=0; T z=0;
			int i = 0; T rho = 0;
			std::string meshFileName, material;
			XMLreader r(this->c->parameterXmlFileName);
			r["obstacle"]["obstacleStartX"].read(x);
			r["obstacle"]["obstacleStartY"].read(y);
			r["obstacle"]["obstacleStartZ"].read(z);
			r["obstacle"]["meshFileName"].read(meshFileName);
			r["obstacle"]["referenceDirection"].read(this->referenceDirection);
			r["obstacle"]["material"].read(material);
			r["obstacle"]["density"].read(rho);
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] MeshFileName ="+meshFileName;
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->position = Point<T>(x,y,z);
			this->velocity[0] = 0;	this->velocity[1] = 0;	this->velocity[2] = 0;
			this->acceleration[0] = 0;	this->acceleration[1] = 0;	this->acceleration[2] = 0;
			this->rotation[0] = 0;	this->rotation[1] = 0; this->rotation[2] = 0;
			this->rotationalVelocity[0] = 0;	this->rotationalVelocity[1] = 0;	this->rotationalVelocity[2] = 0;
			this->rotationalAcceleration[0] = 0;	this->rotationalAcceleration[1] = 0;	this->rotationalAcceleration[2] = 0;
			#ifdef PLB_MPI_PARALLEL
				if(global::mpi().isMainProcessor()){
					this->triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
					global::mpiData().sendTriangleSet<T>(this->triangleSet);
				}
				else{ this->triangleSet = global::mpiData().receiveTriangleSet<T>(); }
			#else
				this->triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
			#endif
			this->flowType = voxelFlag::outside;
			this->density = rho;
			this->volume = this->getVolume();
			this->mass = this->density * this->volume;
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Number of triangles in Mesh = "+std::to_string(this->triangleSet.getTriangles().size());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Done Initializing Obstacle";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}
// Methods
	Obstacle &getCenter(){ // Calculates the center of the obstacle.
		Cuboid<T> cuboid = this->triangleSet.getBoundingCuboid();
		Array<T,3>	lowerLeftCorner = cuboid.lowerLeftCorner;
		Array<T,3>	upperRightCorner = cuboid.upperRightCorner;
		T lowerBound = 0;
		T upperBound = 0;
		lowerBound = std::min(lowerLeftCorner[0],upperRightCorner[0]);
		upperBound = std::max(lowerLeftCorner[0],upperRightCorner[0]);
		this->center.x = (upperBound-lowerBound)/2;
		lowerBound = std::min(lowerLeftCorner[1],upperRightCorner[1]);
		upperBound = std::max(lowerLeftCorner[1],upperRightCorner[1]);
		this->center.y = (upperBound-lowerBound)/2;
		lowerBound = std::min(lowerLeftCorner[2],upperRightCorner[2]);
		upperBound = std::max(lowerLeftCorner[2],upperRightCorner[2]);
		this->center.z = (upperBound-lowerBound)/2;
		return *this;
	}
	T getVolume(){
		T volume = 0;
		this->getCenter();
		std::vector<Array<Array<T,3>,3> > triangles = this->triangleSet.getTriangles();
		for(int i =0; i<triangles.size(); i++){
			Pyramid<T> p(triangles[i],center);
			volume += p.volume();
		}
		return volume;
	}
	// Function to Move the Obstacle through the Fluid
	void move(const plint& dt, std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice){
		Array<T,3> coord = this->mesh->getPhysicalLocation();
		Array<T,3> fluidForce = this->boundaryCondition->getForceOnObject();
		acceleration[0] = fluidForce[0]/this->mass;
		acceleration[1] = fluidForce[1]/this->mass;
		acceleration[2] = fluidForce[2]/this->mass;
		acceleration[2] -= c->gravitationalAcceleration;
		velocity[0] += acceleration[0]*dt;
		velocity[1] += acceleration[1]*dt;
		velocity[2] += acceleration[2]*dt;
		Array<T,3> ds(velocity[0]*dt, velocity[1]*dt, velocity[2]*dt);
		coord[0] += ds[0]; coord[1] += ds[1]; coord[2] += ds[2];
		this->mesh->setPhysicalLocation(coord);
		reparallelize(*lattice);
	}
	// Function to Move Obstacle to it's starting position
	void move(){
		Array<T,3> vec;
		vec[0] = 0;
		vec[1] = 0;
		vec[2] = 0;
		this-> triangleSet.translate(vec);
	}

	void setMesh(DEFscaledMesh<T>* fromMesh){
		this->mesh = fromMesh;
	}

	void setBoundary(OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* fromBoundary){
		this->boundaryCondition = fromBoundary;
	}
// Attributes
	plint referenceDirection;
	int flowType;
	T density, mass, volume, temperature;
	Point<T> center, position;
	Array<T,3> rotation, velocity, rotationalVelocity, acceleration, rotationalAcceleration;
	TriangleSet<T> triangleSet;
	DEFscaledMesh<T>* mesh;
	OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* boundaryCondition;
private:
	bool master;
	const Constants<T,BoundaryType>* c;
};
template<typename T, class BoundaryType>
Obstacle<T,BoundaryType> *Obstacle<T,BoundaryType>::o=0;

template<typename T, class BoundaryType>
class Fluid{
// Default Constructor from Input XML file
public:
	explicit Fluid(Constants<T,BoundaryType>* _c){
		try{
			this->c = _c;
			this->f = this;
			std::string material;
			int i = 0;
			XMLreader r(this->c->parameterXmlFileName);
			r["fluid"]["material"].read(material);
			this->temperature = c->initialTemperature;
			if(material.compare("FLiNaK")==0){i = 1;}
			switch(i){
				case(0):	throw "Fluid Material not Properly Defined! Modify input parameters.xml";
				case(1):	this->density = 2579.3 - 0.6237 * this->temperature;
			}
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}
// Default Destructor
	~Fluid(){delete c;}
// Methods
// Attributes
	T temperature, density;
private:
	Constants<T,BoundaryType>* c;
	const Fluid<T,BoundaryType>* f;
};

template<typename T, class BoundaryType, class SurfaceData>
class Variables{ // Singleton Class
private:
	static Variables<T,BoundaryType,SurfaceData>* v;
	Variables(){ // Default Constructor
		this->resolution = 0; this->gridLevel=0; this->reynolds=0;
		this->p.reset(nullptr);
		this->location = Array<T,3>();
		this->dx = 1;
		this->dt = 1;
		this->iter = 0;
		this->nprocs = 0;
		this->nprocs_side = 0;
		this->first = true;
		this->master = global::mpi().isMainProcessor();
		this->nprocs = global::mpi().getSize();
		this->nprocs_side = (int)cbrt(this->nprocs);
		this->c = Constants<T,BoundaryType>::getInstance();
		this->resolution = this->c->referenceResolution;
		this->w = Wall<T,BoundaryType>::getInstance();
		this->o = Obstacle<T,BoundaryType>::getInstance();
	}
	~Variables(){// Default Destructor
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Variables DESTRUCTOR was called";
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
			throw std::runtime_error("Variables Destructor was Called");
		#endif
		delete c, v, w, o;
	}
public:
// Methods
	static Variables* getInstance(){ if(!v || v == nullptr){ v = new Variables<T,BoundaryType,SurfaceData>();} return v; }

	void update(const plint& _gridLevel, const plint& _reynolds){
		try{
			#ifdef PLB_DEBUG
				std::string mesg ="[DEBUG] Updating Resolution";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->gridLevel = _gridLevel;
			this->resolution = this->c->referenceResolution * util::twoToThePowerPlint(_gridLevel);
			plint scaled_u0lb = this->c->u0lb * util::twoToThePowerPlint(_gridLevel);
			if(c->test){this->reynolds = c->testRe;}
			else{this->reynolds = _reynolds;}
			this->p.reset(new IncomprFlowParam<T>(scaled_u0lb,this->reynolds,this->resolution,1,1,1));
			this->dx = this->p->getDeltaX();
			this->dt = this->p->getDeltaT();
			this->scalingFactor = (T)(this->resolution)/this->dx;
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Resolution="+std::to_string(this->resolution);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Done Updating Resolution";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}

	bool checkConvergence(){
		util::ValueTracer<T> tracer(this->p->getLatticeU(), this->p->getDeltaX(), this->c->epsilon);
		return tracer.hasConverged();
	}

	std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > makeParallel(std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > lattice){
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Parallelizing Lattice ";
			if(this->master){std::cout << mesg << std::endl;}
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
		ParallellizeByCubes3D parallel(domains, box, this->nprocs_side, this->nprocs_side, this->nprocs_side);
		parallel.parallelize();
		#ifdef PLB_DEBUG
			mesg = "[DEBUG] Done Parallelizing Lattice time="+std::to_string(global::timer("parallel").getTime());
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
			global::timer("parallel").stop();
		#endif
		return std::move(lattice);
	}

	std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > getLattice(){
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Lattices";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice = nullptr;
			lattice = getBoundaryCondition(true, std::move(lattice));
			lattice = getBoundaryCondition(false, std::move(lattice));
			lattice->toggleInternalStatistics(false);
			lattice = makeParallel(std::move(lattice));
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Creating Lattices";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			return std::move(lattice);
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}

	std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > getBoundaryCondition(bool wall,
							std::unique_ptr< MultiBlockLattice3D<T,Descriptor> > lattice){
		try{
			plint referenceDirection = 0; int flowType = 0; TriangleSet<T> triangleSet;
			if(wall){
				triangleSet = this->w->triangleSet; referenceDirection = this->w->referenceDirection; flowType = this->w->flowType;
			}
			else{
				triangleSet = this->o->triangleSet; referenceDirection = this->o->referenceDirection; flowType = this->o->flowType;
			}
			#ifdef PLB_DEBUG
				std::string mesg="[DEBUG] Number of Triangles= "+std::to_string(triangleSet.getTriangles().size());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg="[DEBUG] Reference Direction= "+std::to_string(referenceDirection);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg="[DEBUG] FlowType= "+std::to_string(flowType);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg = "[DEBUG] Lattice Addres= "+ adr_string(&lattice);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").start();
			#endif
			DEFscaledMesh<T>* mesh = nullptr;
			mesh = new DEFscaledMesh<T>(triangleSet, this->resolution, referenceDirection, c->margin, c->extraLayer);
			if(wall){ this->w->setMesh(mesh); }else{ this->o->setMesh(mesh); }
			TriangleBoundary3D<T>* triangleBoundary = nullptr;
			triangleBoundary = new TriangleBoundary3D<T>(*mesh,true);
			this->location = triangleBoundary->getPhysicalLocation();
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] TriangleSet address= "+adr_string(&triangleSet);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Mesh address= "+adr_string(mesh);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] TriangleBoundary address= "+adr_string(triangleBoundary);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
			delete mesh;

			triangleBoundary->getMesh().inflate();
			VoxelizedDomain3D<T>* voxelizedDomain = nullptr;
			voxelizedDomain = new VoxelizedDomain3D<T>(*triangleBoundary, flowType, c->extraLayer, c->borderWidth,
															c->envelopeWidth, c->blockSize, this->gridLevel, c->dynamicMesh);
			std::cout << "BUGGG" << std::endl;
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] VoxelizedDomain address= "+adr_string(voxelizedDomain);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
			delete triangleBoundary;

			std::cout << "Creating MultiBlockLattice" << std::endl;
			if(first){ lattice =  generateMultiBlockLattice<T,Descriptor>(voxelizedDomain->getVoxelMatrix(),
												c->envelopeWidth, new IncBGKdynamics<T,Descriptor>(this->p->getOmega()));
			}
			else{*lattice += *generateMultiBlockLattice<T,Descriptor>(voxelizedDomain->getVoxelMatrix(), c->envelopeWidth,
								new IncBGKdynamics<T,Descriptor>(this->p->getOmega()));
			}
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] IncomprFlowParam address= "+adr_string(&p);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Lattice address= "+adr_string(&lattice);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			BoundaryProfiles3D<T,SurfaceData>* profiles = new BoundaryProfiles3D<T,SurfaceData>();
			TriangleFlowShape3D<T,SurfaceData>* flowShape = new TriangleFlowShape3D<T,SurfaceData>(voxelizedDomain->getBoundary(),
				*profiles);
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Profiles address= "+adr_string(profiles);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] FlowShape address= "+adr_string(flowShape);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif
			delete profiles;

			GuoOffLatticeModel3D<T,Descriptor>* model = new GuoOffLatticeModel3D<T,Descriptor>(flowShape, flowType);
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] Model address= "+adr_string(model);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").restart();
			#endif

			OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* boundaryCondition = nullptr;
			boundaryCondition = new OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>(model,*voxelizedDomain, *lattice);
			#ifdef PLB_DEBUG
				mesg ="[DEBUG] BoundaryCondition address= "+adr_string(boundaryCondition);
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				mesg ="[DEBUG] Elapsed Time="+std::to_string(global::timer("boundary").getTime());
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
				global::timer("boundary").stop();
			#endif
			delete flowShape;
			delete model;
			if(wall){ this->w->setBoundary(boundaryCondition);}else{ this->o->setBoundary(boundaryCondition); }
			this->first = false;
			return std::move(lattice);
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}

	std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > saveFields(std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > lattice){
		if(this->iter % this->p->nStep(c->imageSave) == 0){
			lattice->toggleInternalStatistics(true);
			this->boundingBox = lattice->getMultiBlockManagement().getBoundingBox();
			MultiTensorField3D<T,3> v(this->boundingBox.getNx(), this->boundingBox.getNy(), this->boundingBox.getNz());
			computeVelocity(*lattice, v, this->boundingBox);
			this->velocity.push_back(v);
		}
		lattice->toggleInternalStatistics(false);
		return std::move(lattice);
	}

// Attributes
	plint resolution, gridLevel, reynolds, dx, dt, iter;
	Array<T,3> location;
	Box3D boundingBox;
	double time, scalingFactor;
	std::vector<MultiTensorField3D<double,3> > velocity;
private:
	int nprocs, nprocs_side;
	bool master, first;
	const Constants<T,BoundaryType>* c;
	Obstacle<T,BoundaryType>* o;
	Wall<T,BoundaryType>* w;
	std::unique_ptr<IncomprFlowParam<T> > p;
};
template<typename T, class BoundaryType, class SurfaceData>
Variables<T,BoundaryType,SurfaceData> *Variables<T,BoundaryType,SurfaceData>::v=0;

template<typename T, class BoundaryType, class SurfaceData>
class Output{
private:
	static Output<T,BoundaryType,SurfaceData>* out;
	Output(){
		this->c = Constants<T,BoundaryType>::getInstance();
		this->v = Variables<T,BoundaryType,SurfaceData>::getInstance();
		this->master = global::mpi().isMainProcessor();
		bool parallel = false;
		#ifdef PLB_MPI_PARALLEL
			parallel = true;
		#endif
		global::log().init(parallel);
	}
	~Output(){
		delete this->c;
		delete this->v;
		delete this->out;
	}
public:
	// Methods
	static Output* getInstance(){if(!out || out == nullptr){ out = new Output<T,BoundaryType,SurfaceData>();} return out; }

	void elapsedTime(){if(c->test){std::thread(timeLoop);}}

	void timeLoop(){
		double elapsed = 0;
		double testTime = c->testTime*60;
		bool running = true;
		while(running){
			elapsed = this->timer.getTime() - this->startTime;
			if(elapsed > testTime){ throw std::runtime_error("Test Time Limit Exceeded!"); running = false; }
			std::this_thread::sleep_for(std::chrono::minutes(1));
		}
	}

	void writeGifs(MultiBlockLattice3D<T,Descriptor>& lattice, plint iter)
	{
		const plint imSize = 600;
		const plint nx = lattice.getNx();
		const plint ny = lattice.getNy();
		const plint nz = lattice.getNz();
		Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
		ImageWriter<T> imageWriter("leeloo");
		imageWriter.writeScaledGif(createFileName("u", iter, 6),*computeVelocityNorm(lattice, slice),imSize, imSize );
	}
	void writeImages(){
		/*T dx = this->c->parameters[this->v->gridLevel][this->v->reynolds]->getDeltaX();
		T dt = this->c->parameters[this->v->gridLevel][this->v->reynolds]->getDeltaT();
		int n = v->velocity.size();
		for(int i=0; i<n; i++){
			VtkImageOutput3D<T> vtkOut(createFileName("vtk",i,n), dx);
			MultiTensorField3D<T,3> v_field = v->velocity[i];
			vtkOut.writeData<3,float>(v_field, name, dx/dt);
		}*/
	}
	void startMessage(){
		time_t rawtime;
		struct tm * timeinfo;
		time (&rawtime);
		timeinfo = localtime (&rawtime);
		if(this->master){std::cout<<"SIMULATION START "<< asctime(timeinfo) << std::endl;}
		#ifdef PLB_MPI_PARALLEL
		if(this->master){
			int size = plb::global::mpi().getSize();
			std::cout<<"NUMBER OF MPI PROCESSORS="<< size << std::endl;
			//if((int)cbrt(size) % 1){ throw std::runtime_error("Number of MPI Processess must satisfy Cubic Root");}
			std::string imaster =  this->master ? " YES " : " NO ";
			std::cout<<"Is this the main process?"<< imaster << std::endl;
		}
		#endif
		std::string outputDir = "./tmp/";
		global::directories().setOutputDir(outputDir);// Set output DIR w.r.t. to current DIR
	}
	void simMessage(){
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Creating Timer";
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		this->timer.start();		// Start Timer
		this->startTime = this->timer.getTime();
		#ifdef PLB_DEBUG
			mesg = "[DEBUG] Timer Started";
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
			if(c->test){mesg = "[DEBUG] Starting Test";}
			else{mesg = "[DEBUG] Starting Normal Run";}
			if(this->master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
	}

	void stopMessage(){
		std::string mesg = "SIMULATION COMPLETE";
		if(this->master){std::cout << mesg << std::endl;}
		global::log(mesg);
		elapsedTime();
		this->timer.stop();
	}
public:

private:
	bool master;
	global::PlbTimer timer;
	double startTime;
	double endTime;
	const Constants<T,BoundaryType>* c;
	const Variables<T,BoundaryType,SurfaceData>* v;
};
template<typename T, class BoundaryType, class SurfaceData>
Output<T,BoundaryType,SurfaceData> *Output<T,BoundaryType,SurfaceData>::out=0;

}// namespace plb

typedef double T;
typedef plb::Array<T,3> BoundaryType;
typedef plb::Array<T,3> SurfaceData;

int main(int argc, char* argv[]) {
	try{
		plb::plbInit(&argc, &argv, true); // Initialize Palabos
		bool master = plb::global::mpi().isMainProcessor();
		plb::Output<T, BoundaryType, SurfaceData>* out = plb::Output<T, BoundaryType, SurfaceData>::getInstance();
		out->startMessage();
		std::string fileName = "";
		plb::global::argv(argc-1).read(fileName);
		plb::Constants<T,BoundaryType>* c = plb::Constants<T,BoundaryType>::getInstance();
		c->initialize(fileName);// Initialize Constants
		plb::Variables<T, BoundaryType, SurfaceData>* v = plb::Variables<T, BoundaryType, SurfaceData>::getInstance();//Initialize Variables
		plb::Wall<T,BoundaryType>* w = plb::Wall<T,BoundaryType>::getInstance();
		w->initialize(); // Initialize Wall
		plb::Obstacle<T,BoundaryType>* o = plb::Obstacle<T,BoundaryType>::getInstance();
		o->initialize(); // Initialize Obstacle
		out->elapsedTime(); // Initialize Test Timer
		for(plb::plint gridLevel = 0; gridLevel<= c->maxGridLevel; gridLevel++){
			for(plb::plint reynolds = c->minRe; reynolds <= c->maxRe; reynolds++){
				v->update(gridLevel,reynolds);
				std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > lattice = v->getLattice();
				bool converged = false;
				for(int i=0; converged == false; i++)
				{
					v->iter++;
					lattice->collideAndStream();
					lattice = v->saveFields(std::move(lattice));
					o->move(v->dt, std::move(lattice));
					//if(out->elapsedTime(){ break; }
					if(v->checkConvergence()){ converged = true; break; }
					if(c->test){ if(v->iter > c->testIter){ break; }}
				}
			}
			#ifdef PLB_DEBUG
				std::string mesg="N collisions="+std::to_string(v->iter);
				if(master){std::cout << mesg << std::endl;}
				plb::global::log(mesg);
				if(master){std::cout<<"Grid Level="+std::to_string(gridLevel);}
				if(master){std::cout << mesg << std::endl;}
				plb::global::log(mesg);
			#endif
		}
		out->writeImages();
		out->stopMessage();
		return 0;																	// Return Process Completed
	}
	catch(const std::exception& e){
		plb::printTrace();															// Call Functions for Full stack trace
		std::string ex = e.what();
		std::string line = std::to_string(__LINE__);
		plb::global::log("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
		return -1;																// Return Error code
	}
};
