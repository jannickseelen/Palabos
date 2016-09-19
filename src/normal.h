#ifndef NORMAL_H
#define NORMAL_H

namespace plb{

template<typename T>
class SurfaceNormal{
public:
    SurfaceNormal(const std::vector<Array<T,3> >& normals_)
        :normals(normals_)
    { }

    Array<T,3> operator()(pluint id)
    {
        return normals[id];
    }

	void update(const ConnectedTriangleSet<T>& triangleSet){
		plint numVertices = triangleSet.getNumVertices();
		normals.clear();
		for (plint iVertex = 0; iVertex < numVertices; iVertex++) {
			T area = 0;
			Array<T,3> unitNormal = Array<T,3>(0,0,0);
			triangleSet.computeVertexAreaAndUnitNormal(iVertex, area, unitNormal);
			normals[iVertex] = unitNormal;
		}
	}
private:
    std::vector<Array<T,3> > normals;
};

}

#endif
