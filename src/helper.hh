#ifndef HELPER_HH
#define HELPER_HH

#include <palabos3D.hh>
#include "myheaders3D.hh"

namespace plb{


	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Helper<T,BoundaryType,SurfaceData,Descriptor>::Helper(){
		master = global::mpi().isMainProcessor();
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Constructing Helper";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		if(objCount == 0){
			if(Constants<T>::objCount == 0)
			{
				constants = Constants<T>();
			}
			if(Wall<T,BoundaryType,SurfaceData,Descriptor>::objCount == 0)
			{
				wall = Wall<T,BoundaryType,SurfaceData,Descriptor>();
			}
			if(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::objCount == 0)
			{
				obstacle = Obstacle<T,BoundaryType,SurfaceData,Descriptor>();
			}
			if(Variables<T,BoundaryType,SurfaceData,Descriptor>::objCount == 0)
			{
				variables = Variables<T,BoundaryType,SurfaceData,Descriptor>();
			}
			if(Output<T,BoundaryType,SurfaceData,Descriptor>::objCount == 0)
			{
				output = Output<T,BoundaryType,SurfaceData,Descriptor>();
			}
		}
		objCount++;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Helper<T,BoundaryType,SurfaceData,Descriptor>::~Helper(){
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Helper";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		if(Constants<T>::objCount > 0)
		{
			constants.~Constants<T>();
		}
		if(Wall<T,BoundaryType,SurfaceData,Descriptor>::objCount > 0)
		{
			wall.~Wall<T,BoundaryType,SurfaceData,Descriptor>();
		}
		if(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::objCount > 0)
		{
			obstacle.~Obstacle<T,BoundaryType,SurfaceData,Descriptor>();
		}
		if(Variables<T,BoundaryType,SurfaceData,Descriptor>::objCount > 0)
		{
			variables.~Variables<T,BoundaryType,SurfaceData,Descriptor>();
		}
		if(Output<T,BoundaryType,SurfaceData,Descriptor>::objCount > 0)
		{
			output.~Output<T,BoundaryType,SurfaceData,Descriptor>();
		}
		objCount--;
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Helper<T,BoundaryType,SurfaceData,Descriptor>::initialize(const std::string& fileName){
		if(Constants<T>::objCount == 1)
		{
			constants.initialize(fileName);
		}
		if(Wall<T,BoundaryType,SurfaceData,Descriptor>::objCount == 1)
		{
			wall.initialize();
		}
		if(Obstacle<T,BoundaryType,SurfaceData,Descriptor>::objCount == 1)
		{
			obstacle.initialize();
		}
		if(Variables<T,BoundaryType,SurfaceData,Descriptor>::objCount == 1)
		{
			variables.initialize();
		}
		if(Output<T,BoundaryType,SurfaceData,Descriptor>::objCount == 1)
		{
			output.initialize();
		}
	}


} // namespace plb

#endif // HELPER_HH
