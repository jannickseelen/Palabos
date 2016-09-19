#ifndef VELOCITY_H
#define VELOCITY_H

#include <palabos3D.h>
#include "myheaders3D.h"

namespace plb{

template<typename T>
class SurfaceVelocity{
public:
	static int objCount;

	SurfaceVelocity();

	Array<T,3> operator()(pluint id);

	void initialize(const T& mass, const T& g, const T& rho, const T& dt, const T& dx, const ConnectedTriangleSet<T>& triangleSet);

	Array<T,3> getCG(std::vector<Array<T,3> > vertexList);

	Array<T,3> getArm(const Array<T,3>& p1, const Array<T,3>& p2);

	Array<T,6> getMomentOfInertia(const Array<T,3>& cg, const ConnectedTriangleSet<T>& triangleSet);

	Array<T,3> getAlpha(const Array<T,3>& M, const Array<T,6>& I);

	Array<T,3> update(const T& timeLB, const Array<T,3>& force, const Array<T,3>& torque, ConnectedTriangleSet<T>& triangleSet);
// Attributes
private:
	static bool master;
	static bool rotation;
	static std::vector<Array<T,3> > verticesVelocity;
	static std::vector<Array<T,3> > forceList;
	static std::vector<Array<T,3> > torque;
	static std::vector<Array<T,3> > location;
	static std::vector<Array<T,3> > acceleration;
	static std::vector<Array<T,3> > velocity;
	static std::vector<Array<T,3> > angular_acceleration;
	static std::vector<Array<T,3> > angular_velocity;
	static std::vector<T> time;
	static T rho, mass, g, dx, dt;
};

// Initializers
template<typename T>
int SurfaceVelocity<T>::objCount= 0;

template<typename T>
bool SurfaceVelocity<T>::master= false;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::verticesVelocity;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::forceList;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::torque;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::location;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::acceleration;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::velocity;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::angular_acceleration;

template<typename T>
std::vector<Array<T,3> > SurfaceVelocity<T>::angular_velocity;

template<typename T>
std::vector<T> SurfaceVelocity<T>::time;

template<typename T>
T SurfaceVelocity<T>::rho = (T)0;

template<typename T>
T SurfaceVelocity<T>::mass = (T)0;

template<typename T>
T SurfaceVelocity<T>::g = (T)9.81;

template<typename T>
T SurfaceVelocity<T>::dt = (T)1;

template<typename T>
T SurfaceVelocity<T>::dx = (T)1;

}// namespace plb


#endif //VELOCITY_H
