#ifndef MPI_DATA_MANAGER_H
#define MPI_DATA_MANAGER_H

#include "core/globalDefs.h"
#include <parallelism/mpiManager.h>
#include <atomicBlock/dataField3D.h>
#include <offLattice/triangleSet.h>

#ifdef PLB_MPI_PARALLEL
#include <mpi.h>
#include <vector>
#include <string>
#endif

namespace plb{
	namespace global{

	#ifdef PLB_MPI_PARALLEL
	class MpiDataManager{
	public:
		template<typename T>
		void sendScalarField3D(const ScalarField3D<T>& field, const Box3D& fromDomain);

		template<typename T>
		void receiveScalarField3D(ScalarField3D<T>& field, const Box3D& fromDomain,const int& fromId) const;

		template<typename T>
		TriangleSet<T> receiveTriangleSet();

		template<typename T>
		void sendTriangleSet(const TriangleSet<T>& triangles);

		void sendReceiveDomains(const bool& master, std::vector<Box3D>& mpiDomains);

		std::vector<Box3D> splitDomains(const Box3D& domain);

		MpiDataManager& mpiData(){static MpiDataManager instance; return instance;}
	private:
		void checkDomain(int rank, Box3D domain, const int& line);
		Array<int,2> checkIfCrossed(const int& fMin, const int& fMax, const int& n);
		MpiDataManager();
		~MpiDataManager();
	friend MpiDataManager& mpiData();
	};
	#endif

	#ifndef PLB_MPI_PARALLEL
	class MpiDataManager{
	public:
		template<typename T>
		void sendScalarField3D(const ScalarField3D<T>& field, const Box3D& fromDomain){}

		template<typename T>
		void receiveScalarField3D(ScalarField3D<T>& field, const Box3D& fromDomain, const int& fromId) const{}

		template<typename T>
		TriangleSet<T> receiveTriangleSet(){}

		template<typename T>
		void sendTriangleSet(const TriangleSet<T>& triangles){}

		void sendReceiveDomains(const bool& master, std::vector<Box3D>& mpiDomains){}

		void splitDomains(const Box3D& domain){ }

		MpiDataManager& mpiData(){static MpiDataManager instance; return instance;}
	private:
		void checkDomain(int rank, Box3D domain, const int& line){};
		void checkIfCrossed(const int& fMin, const int& fMax, const int& n){ }
		MpiDataManager();
		~MpiDataManager();
	friend MpiDataManager& mpiData();
	};
	#endif

	inline MpiDataManager& mpiData(){static MpiDataManager instance; return instance;}

	}  // namespace global
}  // namespace plb


#endif  // MPI_DATA_MANAGER_H
