#ifndef NORMAL_H
#define NORMAL_H

namespace plb{

template<typename T>
class SurfaceNormal{
public:
	static int objCount;
    SurfaceNormal();

    Array<T,3> operator()(const pluint& id);

	void update(const TriangleBoundary3D<T>* tb);

private:
	static std::vector<Array<T,3> > normals;
	static bool master;
};


// Initializers
template<typename T>
int SurfaceNormal<T>::objCount= 0;

template<typename T>
std::vector<Array<T,3> > SurfaceNormal<T>::normals;

template<typename T>
bool SurfaceNormal<T>::master= false;

}

#endif
