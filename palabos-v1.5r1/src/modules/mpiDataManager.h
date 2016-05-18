#ifndef MPI_DATA_MANAGER_H
#define MPI_DATA_MANAGER_H

#include "core/globalDefs.h"
#include <parallelism/mpiManager.h>

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

		std::vector<Box3D> splitDomains(const Box3D& domain, const plint& minX, const plint& maxX, const plint& minY, const plint& maxY,
			const plint& minZ, const plint& maxZ);

		MpiDataManager& mpiData(){static MpiDataManager instance; return instance;}
	private:
	MpiDataManager();
	~MpiDataManager();
	friend MpiDataManager& mpiData();
	};

	#else
	class MpiDataManager{
	public:
		template<typename T>
		void sendScalarField3D(const ScalarField3D<T>& field, const Box3D& fromDomain){}

		template<typename T>
		void receiveScalarField3D(ScalarField3D<T>& field, const Box3D& fromDomain, const int& fromId) const{}

		std::vector<Box3D> splitDomains(const Box3D& domain, const plint& minX, const plint& maxX, const plint& minY, const plint& maxY,
			const plint& minZ, const plint& maxZ){};
	private:
	MpiDataManager();
	~MpiDataManager();
	friend MpiDataManager& mpiData();
	};
	#endif

	inline MpiDataManager& mpiData(){static MpiDataManager instance; return instance;}

	}  // namespace global
}  // namespace plb


#endif  // MPI_DATA_MANAGER_H
