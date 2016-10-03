#ifndef NORMAL_HH
#define NORMAL_HH

namespace plb{

	template<typename T>
	SurfaceNormal<T>::SurfaceNormal()
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing SurfaceNormal";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Surface Normal already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T>
	Array<T,3> SurfaceNormal<T>::operator()(const pluint& id)
	{
		return normals[id];
	}

	template<typename T>
	void SurfaceNormal<T>::update(const DEFscaledMesh<T>& mesh){
		plint numVertices = mesh.getMesh().getNumVertices();
		normals.clear();
		normals.resize(numVertices);
		normals.reserve(numVertices);
		const bool weightedArea = false;
		for (plint iVertex = 0; iVertex < numVertices; iVertex++){
			normals[iVertex] =mesh.getMesh().computeVertexNormal(iVertex,weightedArea);
		}
	}

} // NAMESPACE PLB

#endif // VELOCITY_HH
