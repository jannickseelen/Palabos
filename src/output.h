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

	Output();

	~Output();
	// Methods
	static void initialize();

	void elapsedTime();

	void timeLoop();

	void writeGifs(MultiBlockLattice3D<T,Descriptor>& lattice, plint iter);

	void writeImages();

	void startMessage();

	void simMessage();

	void stopMessage();

	static std::unique_ptr<Output<T,BoundaryType,SurfaceData,Descriptor> > out;
private:
	static bool master;
	global::PlbTimer timer;
	static double startTime;
	static double endTime;
};

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Output<T,BoundaryType,SurfaceData,Descriptor>::objCount=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Output<T,BoundaryType,SurfaceData,Descriptor>::master=false;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
double Output<T,BoundaryType,SurfaceData,Descriptor>::startTime=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
double Output<T,BoundaryType,SurfaceData,Descriptor>::endTime=0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
std::unique_ptr<Output<T,BoundaryType,SurfaceData,Descriptor> > Output<T,BoundaryType,SurfaceData,Descriptor>::out(nullptr);

} // namespace plb

#endif // OUTPUT_H
