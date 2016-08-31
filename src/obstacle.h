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

	static Obstacle &getCenter();

	static T getVolume();
	// Function to Move the Obstacle through the Fluid
	static void move();
	// Function to Move Obstacle to it's starting position
	static void moveToStart();

// Attributes
	static plint referenceDirection;
	static int flowType;
	static T density, mass, volume, temperature;
	static Point<T> center, position;
	static Array<T,3> rotation, velocity, rotationalVelocity, acceleration, rotationalAcceleration, location;
	static TriangleSet<T> triangleSet;
	static std::unique_ptr<DEFscaledMesh<T> > mesh;
	static std::unique_ptr<TriangleBoundary3D<T> > tb;
	static std::unique_ptr<VoxelizedDomain3D<T> > vd;
	static std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice;
	static std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > bp;
	static std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > fs;
	static std::unique_ptr<ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor> > model;
	static std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > bc;
	static std::unique_ptr<Obstacle<T,BoundaryType,SurfaceData,Descriptor> > o;
private:
	static bool master;
};

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Obstacle<T,BoundaryType,SurfaceData,Descriptor>::objCount= 0;

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
Point<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::position = Point<T>(0,0,0);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Point<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::center = Point<T>(0,0,0);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
TriangleSet<T> Obstacle<T,BoundaryType,SurfaceData,Descriptor>::triangleSet(LDBL);

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
std::unique_ptr<ExtrapolatedGeneralizedOffLatticeModel3D<T,Descriptor> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::model(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::bc(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<Obstacle<T,BoundaryType,SurfaceData,Descriptor> > Obstacle<T,BoundaryType,SurfaceData,Descriptor>::o(nullptr);

}// namespace plb

#endif //OBSTACLE_H
