#ifndef WALL_H
#define WALL_H

#include <palabos3D.h>
#include "myheaders3D.h"
#include <memory>

namespace plb{

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
class Wall{
private:
public:
	static int objCount;

	Wall();

	~Wall();
// Methods
	static void initialize();

	static Box3D getDomain();

// Attributes
	static bool dynamicMesh;
	static plint referenceDirection;
	static int flowType;
	static T temperature, density;
	static ConnectedTriangleSet<T> triangleSet;
	static Array<T,3> location, center;
	static Box3D domain;
	static std::unique_ptr<DEFscaledMesh<T> > mesh;
	static std::unique_ptr<TriangleBoundary3D<T> > tb;
	static std::unique_ptr<VoxelizedDomain3D<T> > vd;
	static std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice;
	static std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > bp;
	static std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > fs;
	static std::unique_ptr<GuoOffLatticeModel3D<T,Descriptor> > model;
	static std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > bc;
	static std::unique_ptr<Wall<T,BoundaryType,SurfaceData,Descriptor> > w;
private:
	static bool master;
};

// Initializers
template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Wall<T,BoundaryType,SurfaceData,Descriptor>::objCount= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Wall<T,BoundaryType,SurfaceData,Descriptor>::dynamicMesh= false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Wall<T,BoundaryType,SurfaceData,Descriptor>::flowType= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Wall<T,BoundaryType,SurfaceData,Descriptor>::referenceDirection= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Wall<T,BoundaryType,SurfaceData,Descriptor>::temperature= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
T Wall<T,BoundaryType,SurfaceData,Descriptor>::density= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
ConnectedTriangleSet<T> Wall<T,BoundaryType,SurfaceData,Descriptor>::triangleSet = ConnectedTriangleSet<T>(TriangleSet<T>(FLT));

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Wall<T,BoundaryType,SurfaceData,Descriptor>::location = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Wall<T,BoundaryType,SurfaceData,Descriptor>::center = Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Box3D Wall<T,BoundaryType,SurfaceData,Descriptor>::domain = Box3D(0,0,0,0,0,0);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<DEFscaledMesh<T> > Wall<T,BoundaryType,SurfaceData,Descriptor>::mesh(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<TriangleBoundary3D<T> > Wall<T,BoundaryType,SurfaceData,Descriptor>::tb(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<VoxelizedDomain3D<T> > Wall<T,BoundaryType,SurfaceData,Descriptor>::vd(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > Wall<T,BoundaryType,SurfaceData,Descriptor>::lattice(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > Wall<T,BoundaryType,SurfaceData,Descriptor>::bp(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > Wall<T,BoundaryType,SurfaceData,Descriptor>::fs(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<GuoOffLatticeModel3D<T,Descriptor> > Wall<T,BoundaryType,SurfaceData,Descriptor>::model(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > Wall<T,BoundaryType,SurfaceData,Descriptor>::bc(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<Wall<T,BoundaryType,SurfaceData,Descriptor> > Wall<T,BoundaryType,SurfaceData,Descriptor>::w(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Wall<T,BoundaryType,SurfaceData,Descriptor>::master= false;
}// namespace plb


#endif //WALL_H
