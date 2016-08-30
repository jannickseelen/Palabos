#ifndef VARIABLES_H
#define VARIABLES_H

#include <palabos3D.h>
#include "myheaders3D.h"

namespace plb{

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
class Variables{
private:
public:
	static int objCount;

	Variables();

	~Variables();
// Methods
	static void initialize();

	void update(const plint& _gridLevel, const plint& _reynolds);

	bool checkConvergence();

	std::unique_ptr<DEFscaledMesh<T> > createMesh(const TriangleSet<T>& triangleSet, const plint& referenceDirection, const int& flowType);

	std::unique_ptr<TriangleBoundary3D<T> > createTB(const DEFscaledMesh<T>& mesh);

	std::unique_ptr<VoxelizedDomain3D<T> > createVoxels(const TriangleBoundary3D<T>& tb, const int& flowType);

	std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > createLattice(VoxelizedDomain3D<T>& voxelizedDomain);

	std::unique_ptr<BoundaryProfiles3D<T,SurfaceData> > createBP();

	std::unique_ptr<TriangleFlowShape3D<T,SurfaceData> > createFS(const VoxelizedDomain3D<T>& vozelizedDomain,
		const BoundaryProfiles3D<T,SurfaceData>& profile);

	std::unique_ptr<OffLatticeModel3D<T,SurfaceData> > createModel(TriangleFlowShape3D<T,SurfaceData>* flowShape,
		const int& flowType);

	std::unique_ptr<OffLatticeBoundaryCondition3D<T,Descriptor,BoundaryType> > createBC(OffLatticeModel3D<T,SurfaceData>* model,
		VoxelizedDomain3D<T>& vozelizedDomain, MultiBlockLattice3D<T,Descriptor>& lt);

	void join();

	void makeParallel();

	void setLattice();

	void saveFields();

// Attributes
	static plint resolution, gridLevel, reynolds, dx, dt, iter;
	static Array<T,3> location;
	static Box3D boundingBox;
	static double time, scalingFactor;
	static std::vector<MultiTensorField3D<double,3> > velocity;
	static IncomprFlowParam<T> p;
	static std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > lattice;
	static std::unique_ptr<Variables<T,BoundaryType,SurfaceData,Descriptor> > v;
private:
	static int nprocs, nprocs_side;
	static bool master;
	static plint scaled_u0lb;
};

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Variables<T,BoundaryType,SurfaceData,Descriptor>::objCount= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::resolution= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::gridLevel= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::reynolds= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::dx= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::dt= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::iter= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Array<T,3> Variables<T,BoundaryType,SurfaceData,Descriptor>::location= Array<T,3>();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
Box3D Variables<T,BoundaryType,SurfaceData,Descriptor>::boundingBox= Box3D();

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
double Variables<T,BoundaryType,SurfaceData,Descriptor>::time= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
double Variables<T,BoundaryType,SurfaceData,Descriptor>::scalingFactor= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Variables<T,BoundaryType,SurfaceData,Descriptor>::nprocs= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Variables<T,BoundaryType,SurfaceData,Descriptor>::nprocs_side= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Variables<T,BoundaryType,SurfaceData,Descriptor>::scaled_u0lb= 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Variables<T,BoundaryType,SurfaceData,Descriptor>::master= false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
IncomprFlowParam<T> Variables<T,BoundaryType,SurfaceData,Descriptor>::p = IncomprFlowParam<T>(scaled_u0lb,reynolds,resolution,1,1,1);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<MultiBlockLattice3D<T,Descriptor> > Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::vector<MultiTensorField3D<double,3> > Variables<T,BoundaryType,SurfaceData,Descriptor>::velocity;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<Variables<T,BoundaryType,SurfaceData,Descriptor> >	Variables<T,BoundaryType,SurfaceData,Descriptor>::v(nullptr);

} // namespace plb

#endif // VARIABLES_H
