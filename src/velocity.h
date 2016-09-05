#ifndef VELOCITY_H
#define VELOCITY_H

#include <palabos3D.h>
#include "myheaders3D.h"

namespace plb{

template<typename T>
class SurfaceVelocity{
private:
	static int objCount;
public:
	SurfaceVelocity();

	Array<T,3> operator()(Array<T,3> const& pos){ return velocity.back();}

	void initialize(const Array<T,3>& start, const T& mass, const T& g);

	Array<T,3> update(const T& timeLB, Array<T,3> force);
// Attributes
private:
	static bool master;
	static std::vector<Array<T,3> > location;
	static std::vector<Array<T,3> > force;
	static std::vector<Array<T,3> > acceleration;
	static std::vector<Array<T,3> > velocity;
	static std::vector<T> time;
	static T mass;
	static T g;
};

// Initializers
template<typename T>
int SurfaceVelocity<T>::objCount= 0;

template<typename T>
bool SurfaceVelocity<T>::master= false;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::location;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::force;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::acceleration;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::velocity;

template<typename T>
std::vector<T> SurfaceVelocity<T>::time;

template<typename T>
T SurfaceVelocity<T>::mass = (T)0;

template<typename T>
T SurfaceVelocity<T>::g = (T)0;

}// namespace plb


#endif //VELOCITY_H
