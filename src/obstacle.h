#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <palabos3D.h>
#include "myheaders3D.h"
#include <memory>

namespace plb{

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
class Obstacle{
private:
public:
	static int objCount;

	Obstacle();
	// Default Destructor

	~Obstacle();
// Methods
	void initialize();

	static void getCenter();

	static void getVolume();

	// Function to Move Obstacle to it's starting position
	static void moveToStart();

	static Array<T,3> getForce();

	// Function to Move the Obstacle through the Fluid
	static void move();

// Attributes
	static SurfaceVelocity<T> surfaceVelocity;
	static bool dynamicMesh, firstMove;
	static plint referenceDirection;
	static int flowType;
	static T density, mass, volume, temperature, numTriangles, numVertices, g;
	static std::vector<Array<T,3> > vertices;
	static std::vector<Array<T,3> > unitNormals;
    static std::vector<T> areas;
	static Point<T> center, position;
	static Array<T,3> rotation, velocity, rotationalVelocity, acceleration, rotationalAcceleration, location;
	static ConnectedTriangleSet<T> triangleSet;
	static std::unique_ptr<DEFscaledMesh<T> > mesh;
	static std::unique_ptr<TriangleBoundary3D<T> > tb;
	static std::unique_ptr<VoxelizedDomain3D<T> > vd;
	static std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice;
	static std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > bp;
	static std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > fs;
	static std::unique_ptr<GuoOffLatticeModel3D<T,Descriptor> > model;
	static std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > bc;
	static std::unique_ptr<SurfaceVelocity<T> > velocityFunc;
	static std::unique_ptr<SurfaceNormal<T> > normalFunc;
	static std::unique_ptr<Obstacle<T,BoundaryType,SurfaceData,Descriptor> > o;
private:
	static Array<T,3> rotation_LB, velocity_LB, rotationalVelocity_LB, acceleration_LB, rotationalAcceleration_LB, location_LB;
	static bool master;
};

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Obstacle<T,BoundaryType,SurfaceData,Descriptor>::objCount= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
SurfaceVelocity<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::surfaceVelocity;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Obstacle<T,BoundaryType,SurfaceData,Descriptor>::dynamicMesh= false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Obstacle<T,BoundaryType,SurfaceData,Descriptor>::firstMove= true;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Obstacle<T,BoundaryType,SurfaceData,Descriptor>::referenceDirection= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Obstacle<T,BoundaryType,SurfaceData,Descriptor>::flowType= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::mass= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::volume= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::density= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::temperature= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::numTriangles= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::numVertices= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,SurfaceData,Descriptor>::g= 9.81;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::vector<Array<T,3> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vertices;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::vector<Array<T,3> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::unitNormals;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::vector<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::areas;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::acceleration = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::rotation = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::rotationalVelocity = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::rotationalAcceleration = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::velocity = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::location = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Obstacle<T,BoundaryType,SurfaceData,Descriptor>::master = false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::acceleration_LB = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::rotation_LB = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::rotationalVelocity_LB = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::rotationalAcceleration_LB = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::velocity_LB = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::location_LB = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Point<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::position = Point<T>(0,0,0);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Point<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::center = Point<T>(0,0,0);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
ConnectedTriangleSet<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::triangleSet = ConnectedTriangleSet<T>(TriangleSet<T>(FLT));

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<DEFscaledMesh<T> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::mesh(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<TriangleBoundary3D<T> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::tb(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<VoxelizedDomain3D<T> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::vd(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::lattice(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bp(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::fs(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<GuoOffLatticeModel3D<T,Descriptor> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::model(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bc(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<SurfaceVelocity<T> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::velocityFunc(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<SurfaceNormal<T> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::normalFunc(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<Obstacle<T,BoundaryType,SurfaceData,Descriptor> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::o(nullptr);

}// namespace plb

#endif //OBSTACLE_H
