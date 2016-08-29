#ifndef GEO_HH
#define GEO_HH

#include "geo.h"

namespace plb{

/* ####################################################################################################################################
 * 						CLASS POINT
 * ####################################################################################################################################*/

	template<typename T>
	Point<T>::Point(){
		this->x = 0;
		this->y = 0;
		this->z = 0;
	}

	template<typename T>
	Point<T>::Point(const T &_x, const T &_y, const T &_z):x(_x),y(_y),z(_z){}

	template<typename T>
	Point<T>::Point(const Array<T,3> &array):x(array[0]),y(array[1]),z(array[2]){}

	template<typename T>
	Point<T>::Point(const Point &p):x(p.x),y(p.y),z(p.z){}

	template<typename T>
	T Point<T>::distance(const Point &p){
		T d = 0;
		d = sqrt(pow((this->x - p.x),2) + pow((this->y - p.y),2) + pow((this->z - p.z),2));
		return d;
	}

/* ####################################################################################################################################
 * 						CLASS TRIANGLE
 * ####################################################################################################################################*/

	template<typename T>
	Triangle<T>::Triangle(const Point<T> &_a, const Point<T> &_b, const Point<T> &_c):a(_a),b(_b),c(_c){}

	template<typename T>
	Triangle<T>::Triangle(const Array<Array<T,3>,3> &array){
			this->a = Point<T>(array[0]);
			this->b = Point<T>(array[1]);
			this->c = Point<T>(array[2]);
	}

	template<typename T>
	Point<T> Triangle<T>::getCentroid(){
		Point<T> center(1/3*(this->a.x + this->b.x + this->c.x),
		1/3*(this->a.y + this->b.y + this->c.y), 1/3*(this->a.z + this->b.z + this->c.z));
		return center;
	}

	template<typename T>
	double Triangle<T>::area() { // Heron's formula
		double area = 0; double AB=0; double BC=0; double CA=0; double s = 0;
		AB = this->a.distance(this->b);
		BC = this->b.distance(this->c);
		CA = this->c.distance(this->a);
		s = (AB + BC + CA)/2;
		area = sqrt(s*(s-AB)*(s-BC)*(s-CA));
		return area;
	}

/* ####################################################################################################################################
 * 						CLASS PYRAMID
 * ####################################################################################################################################*/

	template<typename T>
	Pyramid<T>::Pyramid(const Triangle<T> &b, const Point<T> &t):base(b),top(t){}

	template<typename T>
	T Pyramid<T>::volume() {
		T volume = 0; T area = 0; T height = 0;
		Point<T> centroid = base.getCentroid();
		area = base.area();
		height = centroid.distance(this->top);
		volume = (1/3)*area*height;
		return volume;
	}

}// namespace plb

#endif // GEO_HH
