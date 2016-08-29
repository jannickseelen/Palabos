#ifndef OBSTACLE_H
#define OBSTACLE_H

#include <palabos3D.h>
#include "myheaders3D.h"
#include <memory>

namespace plb{

template<typename T, class BoundaryType, template<class U> class Descriptor>
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
	static void move(const plint& dt, std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice);
	// Function to Move Obstacle to it's starting position
	static void move();

// Attributes
	static plint referenceDirection;
	static int flowType;
	static T density, mass, volume, temperature;
	static Point<T> center, position;
	static Array<T,3> rotation, velocity, rotationalVelocity, acceleration, rotationalAcceleration, location;
	static TriangleSet<T> triangleSet;
	static std::unique_ptr<TriangleBoundary3D<T> > tb;
	static std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice;
	static std::unique_ptr<DEFscaledMesh<T> > mesh;
	static std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > bc;
	static std::unique_ptr<Obstacle<T,BoundaryType,Descriptor> > o;
private:
	static bool master;
};

template<typename T, class BoundaryType, template<class U> class Descriptor>
int Obstacle<T,BoundaryType,Descriptor>::objCount= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
plint Obstacle<T,BoundaryType,Descriptor>::referenceDirection= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
int Obstacle<T,BoundaryType,Descriptor>::flowType= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,Descriptor>::mass= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,Descriptor>::volume= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,Descriptor>::density= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
T Obstacle<T,BoundaryType,Descriptor>::temperature= 0;

template<typename T, class BoundaryType, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,Descriptor>::acceleration = Array<T,3>();

template<typename T, class BoundaryType, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,Descriptor>::rotation = Array<T,3>();

template<typename T, class BoundaryType, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,Descriptor>::rotationalVelocity = Array<T,3>();

template<typename T, class BoundaryType, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,Descriptor>::rotationalAcceleration = Array<T,3>();

template<typename T, class BoundaryType, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,Descriptor>::velocity = Array<T,3>();

template<typename T, class BoundaryType, template<class U> class Descriptor>
Array<T,3> Obstacle<T,BoundaryType,Descriptor>::location = Array<T,3>();

template<typename T, class BoundaryType, template<class U> class Descriptor>
bool Obstacle<T,BoundaryType,Descriptor>::master = false;

template<typename T, class BoundaryType, template<class U> class Descriptor>
Point<T> Obstacle<T,BoundaryType,Descriptor>::position = Point<T>(0,0,0);

template<typename T, class BoundaryType, template<class U> class Descriptor>
Point<T> Obstacle<T,BoundaryType,Descriptor>::center = Point<T>(0,0,0);

template<typename T, class BoundaryType, template<class U> class Descriptor>
TriangleSet<T> Obstacle<T,BoundaryType,Descriptor>::triangleSet(LDBL);

template<typename T, class BoundaryType, template<class U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > Obstacle<T,BoundaryType,Descriptor>::lattice(nullptr);

template<typename T, class BoundaryType, template<class U> class Descriptor>
std::unique_ptr<DEFscaledMesh<T> > Obstacle<T,BoundaryType,Descriptor>::mesh(nullptr);

template<typename T, class BoundaryType, template<class U> class Descriptor>
std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > Obstacle<T,BoundaryType,Descriptor>::bc(nullptr);

template<typename T, class BoundaryType, template<class U> class Descriptor>
std::unique_ptr<TriangleBoundary3D<T> > Obstacle<T,BoundaryType,Descriptor>::tb(nullptr);

template<typename T, class BoundaryType, template<class U> class Descriptor>
std::unique_ptr<Obstacle<T,BoundaryType,Descriptor> > Obstacle<T,BoundaryType,Descriptor>::o(nullptr);

}// namespace plb

#endif //OBSTACLE_H
