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
#include <exception>
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
#include <atomicBlock/dataProcessingFunctional3D.hh>
#include <core/blockLatticeBase3D.hh>
#include <core/units.h>
#include <core/plbInit.hh>
#include <core/dynamicsIdentifiers.hh>
#include <libraryInterfaces/TINYXML_xmlIO.hh>
#include <io/vtkDataOutput.hh>
#include <io/imageWriter.hh>
#include <algorithm/benchmarkUtil.hh>
#include <modules/mpiDataManager.hh>
// LOCAL INCLUDES
#include "LBMexceptions.hh"
#include "backtrace.hh"

#define Descriptor plb::descriptors::D3Q27Descriptor


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
	Constants(){ //Private Constructor for Static Instance
		this->parameterXmlFileName = ""; this->u0lb=0; this->epsilon=0; this->maxRe=0;this->maxGridLevel=0; this->referenceResolution=0;
		this->margin=0; this->borderWidth=0; this->extraLayer=0; this->blockSize=0; this->envelopeWidth=0; this->initialTemperature=0;
		this->gravitationalAcceleration=0; this->dynamicMesh = false; this->parameters=nullptr; this->master = false;
		this->test = false; this->testRe = 0; this->testTime = 0; this->maxT = 0; this->imageSave = 0;
	}
public:
	// Default Destructor
	~Constants(){
		#ifdef PLB_DEBUG
			if(this->master){std::cout << "[DEBUG] Constants DESTRUCTOR was called" << std::endl;}
		#endif
		if(this->master){throw std::runtime_error("Constants Destructor was Called");}
		delete c;
		delete parameters;
	}
// Methods
	static Constants& getInstance(){static Constants<T,BoundaryType> instance; return instance; }
	Constants& initialize(const std::string& fileName, const bool _master){
		try{
			#ifdef PLB_DEBUG
				if(_master){std::cout << "[DEBUG] Creating Constants" << std::endl;}
			#endif
			this->master = _master;
			this->c = this;
			this->parameterXmlFileName = fileName;
			XMLreader r(fileName);
			r["lbm"]["u0lb"].read(this->u0lb);
			double x = 1;
			double y = 3;
			double cs = sqrt(x / y);
			double mach = this->u0lb / cs;
			double maxMach = 0.1;
			if(mach > maxMach){throw machEx;}
			r["lbm"]["maxRe"].read(this->maxRe);
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
			r["simulation"]["test"].read(this->test);
			r["simulation"]["testRe"].read(this->testRe);
			r["simulation"]["testTime"].read(this->testTime);
			r["simulation"]["maxT"].read(this->maxT);
			r["simulation"]["imageSave"].read(this->imageSave);
			// Initialize a Dynamic 2D array of the flow parameters and following BGKdynamics
			plint row = this->maxGridLevel+1; plint col = 0;
			if(test){ col = 1; }else{ col = this->maxRe+1; }
			this->parameters = new IncomprFlowParam<T>**[row];
			for(plint grid = 0; grid <= this->maxGridLevel; grid++){
				this->parameters[grid] = new IncomprFlowParam<T>*[testRe];
			}
			double ratio = 3;
			// Fill the 2D array with standard values
			for(plint grid = 0; grid <= this->maxGridLevel; grid++){
				T resolution = this->referenceResolution * util::twoToThePowerPlint(grid);
				T scaled_u0lb = this->u0lb * util::twoToThePowerPlint(grid);
				mach = scaled_u0lb / cs;
				if(mach > maxMach){std::cout<<"Local Mach= "<<mach<<"\n"; throw localMachEx;}
				if(resolution == 0){throw resolEx;}
				if(test){
					this->parameters[grid][testRe] = new IncomprFlowParam<T>(scaled_u0lb,testRe,resolution,1,1,1);
					// Check local speed of sound constraint
					T dt = this->parameters[grid][testRe]->getDeltaT();
					T dx = this->parameters[grid][testRe]->getDeltaX();
					if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
				}
				else{
					for(int reynolds = 0; reynolds <= this->maxRe; reynolds++){
						this->parameters[grid][reynolds] = new IncomprFlowParam<T>(scaled_u0lb,reynolds,resolution,1,1,1);
						// Check local speed of sound constraint
						T dt = this->parameters[grid][reynolds]->getDeltaT();
						T dx = this->parameters[grid][reynolds]->getDeltaX();
						if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
					}
				}
			}
			#ifdef PLB_DEBUG
				if(_master){std::cout << "[DEBUG] Done Creating Constants" << std::endl;}
			#endif
			return *this;
		}
		catch(std::exception& e){
			std::cout << "Exception Caught: " << e.what() << "\n";
			throw;
		}
	}
// Methods
// Properties
	std::string parameterXmlFileName;
	T u0lb, epsilon, maxT, imageSave;
	plint testRe, testTime, maxRe, maxGridLevel, referenceResolution, margin, borderWidth, extraLayer, blockSize, envelopeWidth;
	T initialTemperature, gravitationalAcceleration;
	bool dynamicMesh, test;
	IncomprFlowParam<T>*** parameters;
private:
	bool master;
	const Constants<T,BoundaryType>* c;
};


template<typename T, class BoundaryType>
class Wall{
private:
	Wall(){// Private Constructor for static Instance
		this->referenceDirection = 0;
		this->flowType = 0;
		this->temperature = 0; this->density=0;
		this->boundaryCondition = nullptr;
	}
public:
	~Wall(){
		#ifdef PLB_DEBUG
			std::cout << "[DEBUG] Wall DESTRUCTOR was called" << std::endl;
		#endif
		delete c;
		delete boundaryCondition;
		delete w;
	}
// Methods
	static Wall& getInstance(){static Wall<T,BoundaryType> instance;return instance;}
	Wall& initialize(const Constants<T,BoundaryType>* _c, const bool master){
		#ifdef PLB_DEBUG
			if(master){std::cout << "[DEBUG] Initializing Wall" << std::endl;}
		#endif
		this->c = _c;
		this->w = this;
		std::string meshFileName, material;
		int i = 0;
		try{
			XMLreader r(this->c->parameterXmlFileName);
			r["wall"]["meshFileName"].read(meshFileName);
			r["wall"]["referenceDirection"].read(this->referenceDirection);
			r["wall"]["initialTemperature"].read(this->temperature);
			r["wall"]["material"].read(material);
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] MeshFileName =" << meshFileName << std::endl;}
			#endif
		}
		catch (PlbIOException& exception) {
			std::cout << "Error Constructing Wall from: " << this->c->parameterXmlFileName << ": " << exception.what() << std::endl;
			throw;
		}
		this->triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
		this->flowType = voxelFlag::inside;
		this->temperature = this->c->initialTemperature;
		if(material.compare("AL203")==0)
			i = 1;
		switch(i){
			case(0):	throw "Wall Material not Properly Defined! Modify input parameters.xml";
			case(1):	this->density = 3840; //[kg/m3];
		}
		#ifdef PLB_DEBUG
			if(master){std::cout << "[DEBUG] Number of triangles in Mesh = "<< this->triangleSet.getTriangles().size() << std::endl;}
			if(master){std::cout << "[DEBUG] Done Initializing Wall" << std::endl;}
		#endif
		return *this;
	}
// Attributes
	plint referenceDirection;
	int flowType;
	T temperature, density;
	TriangleSet<T> triangleSet;
	OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* boundaryCondition;
private:
	const Constants<T,BoundaryType>* c;
	const Wall<T,BoundaryType>* w;
};

template<typename T, class BoundaryType>
class Obstacle{
private:
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
		this->master = false;
	}
public:
	// Default Destructor
	~Obstacle(){
		#ifdef PLB_DEBUG
			if(master){std::cout << "[DEBUG] Obstacle DESTRUCTOR was called" << std::endl;}
		#endif
		delete c;
		delete boundaryCondition;
		delete o;
	}
	static Obstacle& getInstance(){static Obstacle<T,BoundaryType> instance; return instance; }
	Obstacle& initialize(const Constants<T,BoundaryType>* _c, const bool _master){
		this->master = _master;
		#ifdef PLB_DEBUG
			if(master){std::cout << "[DEBUG] Initializing Obstacle" << std::endl;}
		#endif
		this->c = _c;
		this->o = this;
		T x=0; T y=0; T z=0;
		int i = 0; T rho = 0;
		std::string meshFileName, material;
		try{
			XMLreader r(this->c->parameterXmlFileName);
			r["obstacle"]["obstacleStartX"].read(x);
			r["obstacle"]["obstacleStartY"].read(y);
			r["obstacle"]["obstacleStartZ"].read(z);
			r["obstacle"]["meshFileName"].read(meshFileName);
			r["obstacle"]["referenceDirection"].read(this->referenceDirection);
			r["obstacle"]["material"].read(material);
			r["obstacle"]["density"].read(rho);
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] MeshFileName =" << meshFileName << std::endl;}
			#endif
		}
		catch (PlbIOException& exception) {
			std::cout << "Error Constructing Obstacle from: " << this->c->parameterXmlFileName << ": " << exception.what() << std::endl;
            throw exception;
		}
		this->position = Point<T>(x,y,z);
		this->velocity[0] = 0;	this->velocity[1] = 0;	this->velocity[2] = 0;
		this->acceleration[0] = 0;	this->acceleration[1] = 0;	this->acceleration[2] = 0;
		this->rotation[0] = 0;	this->rotation[1] = 0; this->rotation[2] = 0;
		this->rotationalVelocity[0] = 0;	this->rotationalVelocity[1] = 0;	this->rotationalVelocity[2] = 0;
		this->rotationalAcceleration[0] = 0;	this->rotationalAcceleration[1] = 0;	this->rotationalAcceleration[2] = 0;
		this->triangleSet = TriangleSet<T>(meshFileName, DBL, STL);
		this->flowType = voxelFlag::outside;
		this->density = rho;
		this->volume = this->getVolume();
		this->mass = this->density * this->volume;
		#ifdef PLB_DEBUG
			if(master){std::cout << "[DEBUG] Number of triangles in Mesh = "<< this->triangleSet.getTriangles().size() << std::endl;}
			if(master){std::cout << "[DEBUG] Done Initializing Obstacle" << std::endl;}
		#endif
		return *this;
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
	Obstacle &move(const plint &dt){
		Array<T,3> fluidForce = this->boundaryCondition.getForceOnObject();
		fluidForce[2] -= c->gravitationalAcceleration * this-> mass;
		this-> triangleSet.translate();
		this-> triangleSet.rotate();
		return *this;
	}
	// Function to Move Obstacle to it's starting position
	Obstacle &move(){
		Array<T,3> vec;
		vec[0] = 0;
		vec[1] = 0;
		vec[2] = 0;
		this-> triangleSet.translate(vec);
		return *this;
	}
	Obstacle &getAcceleration(const Array<T,3> &force){}
// Attributes
	plint referenceDirection;
	int flowType;
	T density, mass, volume, temperature;
	Point<T> center, position;
	Array<T,3> rotation, velocity, rotationalVelocity, acceleration, rotationalAcceleration;
	TriangleSet<T> triangleSet;
	OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* boundaryCondition;
private:
	bool master;
	const Obstacle<T,BoundaryType>* o;
	const Constants<T,BoundaryType>* c;
};

template<typename T, class BoundaryType>
class Fluid{
// Default Constructor from Input XML file
public:
	explicit Fluid(Constants<T,BoundaryType>* _c){
		this->c = _c;
		this->f = this;
		std::string material;
		int i = 0;
		try{
			XMLreader r(this->c->parameterXmlFileName);
			r["fluid"]["material"].read(material);
		}
		catch (PlbIOException& exception) {
			std::cout << "Error Constructing Fluid from: " << this->c->parameterXmlFileName << ": " << exception.what() << std::endl;
			throw;
		}
		this->temperature = c->initialTemperature;
		if(material.compare("FLiNaK")==0){i = 1;}
		switch(i){
			case(0):	throw "Fluid Material not Properly Defined! Modify input parameters.xml";
			case(1):	this->density = 2579.3 - 0.6237 * this->temperature;
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
class Variables{
private:
	Variables(){ // Default Constructor
		this->resolution = 0; this->gridLevel=0; this->reynolds=0;
		this->location = Array<T,3>();
		this->first = true; this->master = false;
	}
public:
	~Variables(){// Default Destructor
		#ifdef PLB_DEBUG
			if(master){std::cout << "[DEBUG] Variables DESTRUCTOR was called" << std::endl;}
		#endif
		delete c;
		delete v;
	}
// Methods
	Variables& initialize(const Constants<T,BoundaryType>* _c, const bool _master){
		try{
			this->master = _master;
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Creating Variables" << std::endl;}
			#endif
			this->v = this; // To make sure that there can only be one instance of this class
			this->c = _c;
			this->gridLevel = 0;
			this->reynolds = 0;
			this->resolution = this->c->referenceResolution;
			this->first = true;
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Done Creating Variables" << "\n";}
			#endif
			return *this;
		}
		catch(const std::exception& e){
			std::cout << "Exception Caught in Variable Initializer: " << e.what() << "\n";
			throw;
		}
	}
	static Variables& getInstance(){static Variables<T,BoundaryType,SurfaceData> instance; return instance; }
	void update(const plint& _gridLevel, const plint& _reynolds){
		try{
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Updating Resolution" << std::endl;}
			#endif
			this->gridLevel = _gridLevel;
			this->resolution = this->resolution * util::twoToThePowerPlint(_gridLevel);
			this->dx = this->c->parameters[this->gridLevel][this->reynolds]->getDeltaX();
			this->dt = this->c->parameters[this->gridLevel][this->reynolds]->getDeltaT();
			this->scalingFactor = (T)(this->resolution)/this->dx;
			this->reynolds = _reynolds;
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Resolution=" << this->resolution << std::endl;}
				if(master){std::cout << "[DEBUG] Done Updating Resolution" << std::endl;}
			#endif
		}
		catch(const std::exception& e){
			std::cout << "Exception Caught in update: " << e.what() << "\n";
			throw;
		}
	}
	bool checkConvergence(){
		IncomprFlowParam<T> p = *this->c->parameters[this->gridLevel][this->reynolds];
		util::ValueTracer<T> tracer(p.getLatticeU(), p.getDeltaX(), this->c->epsilon);
		return tracer.hasConverged();
	}
	std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > getLattice(const Wall<T,BoundaryType>& w, const Obstacle<T,BoundaryType>& o){
		try{
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Creating Lattices" << std::endl;}
			#endif
			std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice;
			lattice = getBoundaryCondition(true, w.triangleSet, w.referenceDirection, w.flowType, std::move(lattice));
			lattice = getBoundaryCondition(false, o.triangleSet, o.referenceDirection, o.flowType, std::move(lattice));
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Done Creating Lattices" << std::endl;}
			#endif
			this->boundingBox = lattice->getMultiBlockManagement().getBoundingBox();
			MultiTensorField3D<T,3> v(this->boundingBox.getNx(), this->boundingBox.getNy(), this->boundingBox.getNz());
			computeVelocity(*lattice, v, this->boundingBox);
			this->velocity.push_back(v);
			return lattice;
		}
		catch(const std::exception& e){
			std::cout << "Exception Caught in getLattice: " << e.what() << "\n";
			throw;
		}
	}
	std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > getBoundaryCondition(bool wall, const TriangleSet<T>& triangleSet,
		const plint& referenceDirection, const int& flowType, std::unique_ptr< MultiBlockLattice3D<T,Descriptor> > lattice){
		try{
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Number of Triangles= " << triangleSet.getTriangles().size() <<std::endl;}
				if(master){std::cout << "[DEBUG] Reference Direction= " << referenceDirection <<std::endl;}
				if(master){std::cout << "[DEBUG] FlowType= " << flowType << std::endl;}
				if(master){std::cout << "[DEBUG] Lattice Addres= " << &lattice  << std::endl;}
				if(master){global::timer("boundary").start();}
			#endif
			DEFscaledMesh<T>* mesh = new DEFscaledMesh<T>(triangleSet, this->resolution, referenceDirection, c->margin, c->extraLayer);
			TriangleBoundary3D<T>* triangleBoundary = new TriangleBoundary3D<T>(*mesh,true);
			this->location = triangleBoundary->getPhysicalLocation();
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] TriangleSet address= "<< &triangleSet << std::endl;}
				if(master){std::cout << "[DEBUG] Mesh address= "<< &mesh << std::endl;}
				if(master){std::cout << "[DEBUG] TriangleBoundary address= "<< &triangleBoundary << std::endl;}
				if(master){std::cout << "[DEBUG] Elapsed Time="<<global::timer("boundary").getTime() << std::endl;}
				if(master){global::timer("boundary").restart();}
			#endif
			delete mesh;

			triangleBoundary->getMesh().inflate();
			VoxelizedDomain3D<T>* voxelizedDomain = new VoxelizedDomain3D<T>(*triangleBoundary, flowType, c->extraLayer, c->borderWidth,
															c->envelopeWidth, c->blockSize, this->gridLevel, c->dynamicMesh);
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] VoxelizedDomain address= "<< &voxelizedDomain << std::endl;}
				if(master){std::cout << "[DEBUG] Elapsed Time="<<global::timer("boundary").getTime() << std::endl;}
				if(master){global::timer("boundary").restart();}
			#endif
			delete triangleBoundary;

			IncomprFlowParam<T> p = *this->c->parameters[this->gridLevel][this->reynolds];
			if(first){ lattice =  generateMultiBlockLattice<T,Descriptor>(voxelizedDomain->getVoxelMatrix(),
												c->envelopeWidth, new IncBGKdynamics<T,Descriptor>(p.getOmega()));
			}
			else{*lattice += *generateMultiBlockLattice<T,Descriptor>(voxelizedDomain->getVoxelMatrix(), c->envelopeWidth,
								new IncBGKdynamics<T,Descriptor>(p.getOmega()));
			}
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] IncomprFlowParam address= "<< &p << std::endl;}
				if(master){std::cout << "[DEBUG] Lattice address= "<< &lattice << std::endl;}
				if(master){std::cout << "[DEBUG] Elapsed Time="<< global::timer("boundary").getTime() << std::endl;}
				if(master){global::timer("boundary").restart();}
			#endif

			BoundaryProfiles3D<T,SurfaceData>* profiles = new BoundaryProfiles3D<T,SurfaceData>();
			TriangleFlowShape3D<T,SurfaceData>* flowShape = new TriangleFlowShape3D<T,SurfaceData>(voxelizedDomain->getBoundary(),
				*profiles);
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Profiles address= "<< &profiles << std::endl;}
				if(master){std::cout << "[DEBUG] FlowShape address= "<< &flowShape << std::endl;}
				if(master){std::cout << "[DEBUG] Elapsed Time="<< global::timer("boundary").getTime() << std::endl;}
				if(master){global::timer("boundary").restart();}
			#endif
			delete profiles;

			GuoOffLatticeModel3D<T,Descriptor>* model = new GuoOffLatticeModel3D<T,Descriptor>(flowShape, flowType);
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] Model address= "<< &model << std::endl;}
				if(master){std::cout << "[DEBUG] Elapsed Time="<< global::timer("boundary").getTime() << std::endl;}
				if(master){global::timer("boundary").restart();}
			#endif

			OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>* boundaryCondition = nullptr;
			boundaryCondition = new OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType>(model,*voxelizedDomain, *lattice);
			#ifdef PLB_DEBUG
				if(master){std::cout << "[DEBUG] BoundaryCondition address= "<< &boundaryCondition << std::endl;}
				if(master){std::cout << "[DEBUG] Elapsed Time="<< global::timer("boundary").getTime() << std::endl;}
				if(master){global::timer("boundary").stop();}
			#endif
			delete flowShape;
			delete model;
			if(wall){Wall<T,BoundaryType> w = Wall<T,BoundaryType>::getInstance(); w.boundaryCondition = boundaryCondition; }
			else{Obstacle<T,BoundaryType> o = Obstacle<T,BoundaryType>::getInstance(); o.boundaryCondition = boundaryCondition;}
			this->first = false;
			return lattice;
		}
		catch(const std::exception& e){
			std::cout << "Exception Caught in getBoundaryCondition: " << e.what() << "\n";
			throw;
		}
	}
// Attributes
	plint resolution, gridLevel, reynolds;
	plint dx, dt;
	Array<T,3> location;
	Box3D boundingBox;
	double time, scalingFactor;
	std::vector<MultiTensorField3D<double,3> > velocity;
private:
	bool master;
	bool first;
	const Constants<T,BoundaryType>* c;
	const Variables<T,BoundaryType,SurfaceData>* v;
};

template<typename T, class BoundaryType, class SurfaceData>
class Output{
public:
	Output(const Constants<T,BoundaryType> _c, const Variables<T,BoundaryType,SurfaceData> _v){
		#ifdef PLB_DEBUG
			if(global::mpi().isMainProcessor()){ std::cout << "[DEBUG] Initializing Output" << std::endl; }
		#endif
		this->c = &_c;
		this->v = &_v;
		#ifdef PLB_DEBUG
			if(global::mpi().isMainProcessor()){ std::cout << "[DEBUG] Done Initializing Output" << std::endl; }
		#endif
	}
	// Methods
	bool timeSpent(const plb::global::PlbTimer& timer, const double& startTime){
		double sec = timer.getTime();
		#if PLB_DEBUG
			if(global::mpi().isMainProcessor()){ std::cout << "[DEBUG] Elapsed time in seconds=  " << sec-startTime << std::endl; }
		#endif
		if(c->test){
			if(c->testTime*60 > sec){
				return true;
			}
		}
		return false;
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
private:
	const Constants<T,BoundaryType>* c;
	const Variables<T,BoundaryType,SurfaceData>* v;
};

}// namespace plb

typedef double T;
typedef plb::Array<T,3> BoundaryType;
typedef plb::Array<T,3> SurfaceData;

int main(int argc, char* argv[]) {
	try{
		plb::plbInit(&argc, &argv, true); // Initialize Palabos
		std::string outputDir = "./tmp/";
		bool master = false;
		#ifdef PLB_MPI_PARALLEL
			master = plb::global::mpi().isMainProcessor();
			if(master){
				std::cout<<"SIMULATION START"<< std::endl;
				int size = plb::global::mpi().getSize();
				std::cout<<"NUMBER OF MPI PROCESSORS="<< size << std::endl;
				if((int)cbrt(size) % 2){ throw std::runtime_error("Number of MPI Processess must satisfy Cubic Root");}
				std::string imaster =  master ? " YES " : " NO ";
				std::cout<<"Is this the main process?"<< imaster << std::endl;
				#ifdef PLB_DEBUG
					#define MASTER
				#endif
			}
			// Output the MPI processes whith node identification
			//char processor_name[MPI_MAX_PROCESSOR_NAME];
			//int name_len;
			//MPI_Get_processor_name(processor_name, &name_len);
			//std::cout<<"Processor= "<< processor_name << " ID= "<< plb::global::mpi().getRank() << std::endl;
		#endif
		plb::global::directories().setOutputDir(outputDir);// Set output DIR w.r.t. to current DIR
		std::string fileName = "";
		plb::global::argv(argc-1).read(fileName);
		plb::Constants<T,BoundaryType> c = plb::Constants<T,BoundaryType>::getInstance();
		c.initialize(fileName, master);// Initialize Constants
		plb::Constants<T,BoundaryType>* c_pointer = &c;
		plb::Variables<T, BoundaryType, SurfaceData> v = plb::Variables<T, BoundaryType, SurfaceData>::getInstance();
		v.initialize(c_pointer, master);// Initialize Variables
		plb::Wall<T,BoundaryType> w = plb::Wall<T,BoundaryType>::getInstance();
		w.initialize(c_pointer, master); // Initialize Wall
		plb::Obstacle<T,BoundaryType> o = plb::Obstacle<T,BoundaryType>::getInstance();
		o.initialize(c_pointer, master); // Initialize Obstacle
		#ifdef PLB_DEBUG
			if(master){std::cout<<"[DEBUG] Creating Timer" << std::endl;}
		#endif
		plb::global::PlbTimer timer; timer.start();		// Start Timer
		double startTime = timer.getTime();
		#ifdef PLB_DEBUG
			if(master){std::cout<<"[DEBUG] Timer Started" << std::endl;}
		#endif
		double collisions = 0;
		plb::Output<T, BoundaryType, SurfaceData>* output = new plb::Output<T, BoundaryType, SurfaceData>(c,v); // Initialize Output
		for(plb::plint gridLevel = 0; gridLevel<= c.maxGridLevel; gridLevel++){
			if(c.test == true){
				#ifdef PLB_DEBUG
					if(master){std::cout<<"[DEBUG] Starting Test" << std::endl;}
				#endif
				v.update(gridLevel,c.testRe);								// Update resolution at different gridLevels
				std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > lattice = v.getLattice(w,o);
				bool converged = false;
				for(int i=0; converged == false; i++)
				{
					if(master){if(i % (int)(v.dt * c.imageSave) == 0){ output->writeGifs(*lattice,i);}}
					collisions++;
					lattice->collideAndStream();
					converged = v.checkConvergence();
					if(output->timeSpent(timer, startTime)){ break; }
					if(converged){ break; }
				}
			}
			else{
				#ifdef PLB_DEBUG
					if(master){std::cout<<"[DEBUG] Starting Normal Run" << std::endl;}
				#endif
				for(plb::plint reynolds = 0; reynolds <= c.maxRe; reynolds++){
					v.update(gridLevel,reynolds);								// Update resolution at different gridLevels
					std::unique_ptr<plb::MultiBlockLattice3D<T,Descriptor> > lattice = v.getLattice(w,o);
					bool converged = false;
					for(int i=0; converged == false; i++)
					{
						if(master){if(i % (int)(v.dt * c.imageSave) == 0){ output->writeGifs(*lattice,i);}}
						collisions++;
						lattice->collideAndStream();
						converged = v.checkConvergence();
						if(output->timeSpent(timer, startTime)){ break; }
						if(converged){ break; }
					}
				}
			}
			#ifdef PLB_DEBUG
				if(master){std::cout<<"N collisions=" << collisions << std::endl;}
				if(master){std::cout<<"Grid Level=" << gridLevel << std::endl;}
			#endif
			output->writeImages();
		}
		delete output;
		if(master){std::cout<<"SIMULATION COMPLETE"<< std::endl;}
		if(master){std::cout<< "Total Run Time: " << plb::global::timer("global").getTime() << '\n';}		// Output Elapsed Time
		timer.stop();
		return 0;																	// Return Process Completed
	}
	catch(std::exception& e){
		plb::printTrace();															// Call Functions for Full stack trace
		std::cerr << e.what() << '\n';												// Output exception
		return -1;																// Return Error code
	}
};
