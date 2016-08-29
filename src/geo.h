#ifndef GEO_H
#define GEO_H

#include <core/array.h>

namespace plb{

template<typename T>
class Point{
// Default Constructor
public:
	Point();
	Point(const T &_x, const T &_y, const T &_z);
	Point(const Array<T,3> &array);
	Point(const Point &p);
// Default Destructor
	~Point(){}
// Methods
	T distance(const Point &p);
// Properties
T x;
T y;
T z;
};

template<typename T>
class Triangle{
// Default Constructor
public:
	Triangle(){} // Default constructor that calls the default point constructor for a,b and c
	Triangle(const Point<T> &_a, const Point<T> &_b, const Point<T> &_c);
	Triangle(const Array<Array<T,3>,3> &array);
// Default Destructor
	~Triangle(){};
// Methods
	Point<T> getCentroid();
	double area();
// Properties
	Point<T> a, b, c;
};

template<typename T>
class Pyramid{
// Default Constructor
public:
	Pyramid(){} // Default constructor that calls the default point constructor for a,b and c
	Pyramid(const Triangle<T> &b, const Point<T> &t);
// Default Destructor
	~Pyramid(){};
// Methods
	T volume();
// Properties
	Triangle<T> base;
	Point<T> top;
};
}// Namespace plb

#endif // GEO_H
