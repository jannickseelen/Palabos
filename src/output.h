#ifndef OUTPUT_H
#define OUTPUT_H

#include <palabos3D.h>
#include "myheaders3D.h"

namespace plb{

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
class Output{
private:
public:
	static int objCount;

	explicit Output();

	~Output();
	// Methods
	static void initialize();

	void elapsedTime();

	void timeLoop();

	void writeGif();
private:
	void writeDensity();

	void writeVelocity();

	void writeVorticity();

public:
	void writeImages(const plint& reynolds_, const plint& gridLevel_, const bool& last = false);

	void startMessage();

	void simMessage();

	void stopMessage();

	static std::unique_ptr<Output<T,BoundaryType,SurfaceData,Descriptor> > out;
private:
	static std::unique_ptr<VtkStructuredImageOutput3D<T> > densityOut;
	static std::unique_ptr<VtkStructuredImageOutput3D<T> > velocityOut;
	static std::unique_ptr<VtkStructuredImageOutput3D<T> > vorticityOut;
	static std::unique_ptr<MultiTensorField3D<T,3> > v;
	static std::unique_ptr<MultiTensorField3D<T,3> > w;
	static std::unique_ptr<MultiScalarField3D<T> > r;
	static plint reynolds;
	static plint gridLevel;
	static bool master;
	static bool first;
	static bool last;
	global::PlbTimer timer;
	static double startTime;
	static double endTime;
	static int gifCount;
	static int vtkCount;
};

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Output<T,BoundaryType,SurfaceData,Descriptor>::objCount=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<VtkStructuredImageOutput3D<T> > Output<T,BoundaryType,SurfaceData,Descriptor>::densityOut(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<VtkStructuredImageOutput3D<T> > Output<T,BoundaryType,SurfaceData,Descriptor>::velocityOut(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<VtkStructuredImageOutput3D<T> > Output<T,BoundaryType,SurfaceData,Descriptor>::vorticityOut(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T,3> > Output<T,BoundaryType,SurfaceData,Descriptor>::v(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<MultiTensorField3D<T,3> > Output<T,BoundaryType,SurfaceData,Descriptor>::w(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<MultiScalarField3D<T> > Output<T,BoundaryType,SurfaceData,Descriptor>::r(nullptr);

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Output<T,BoundaryType,SurfaceData,Descriptor>::reynolds=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
plint Output<T,BoundaryType,SurfaceData,Descriptor>::gridLevel=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Output<T,BoundaryType,SurfaceData,Descriptor>::gifCount=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Output<T,BoundaryType,SurfaceData,Descriptor>::vtkCount=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Output<T,BoundaryType,SurfaceData,Descriptor>::master=false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Output<T,BoundaryType,SurfaceData,Descriptor>::first=true;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Output<T,BoundaryType,SurfaceData,Descriptor>::last=false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
double Output<T,BoundaryType,SurfaceData,Descriptor>::startTime=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
double Output<T,BoundaryType,SurfaceData,Descriptor>::endTime=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<Output<T,BoundaryType,SurfaceData,Descriptor> > Output<T,BoundaryType,SurfaceData,Descriptor>::out(nullptr);

} // namespace plb

#endif // OUTPUT_H
