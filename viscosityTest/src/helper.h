#ifndef HELPER_H
#define HELPER_H

#include <Palabos3D.h>
#include "myheaders3D.h"

#include <string>

namespace plb{

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
class Helper{
public:
	static int objCount;
	Helper();
	~Helper();
	void initialize(const std::string& fileName);
private:
//Classes
	Constants<T> constants; //= Constants<T>();
	Wall<T,BoundaryType,SurfaceData,Descriptor> wall; //= Wall<T,BoundaryType,SurfaceData,Descriptor>();
	Obstacle<T,BoundaryType,Descriptor> obstacle; //= Obstacle<T,BoundaryType,Descriptor>();
	Variables<T,BoundaryType,SurfaceData,Descriptor> variables; //= Variables<T,BoundaryType,SurfaceData,Descriptor>();
	Output<T,BoundaryType,SurfaceData,Descriptor> output; //= Output<T,BoundaryType,SurfaceData,Descriptor>();
	static bool master;
};

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
int Helper<T,BoundaryType,SurfaceData,Descriptor>::objCount = 0;

template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
bool Helper<T,BoundaryType,SurfaceData,Descriptor>::master = false;

}//namespace plb

#endif //HELPER_H
