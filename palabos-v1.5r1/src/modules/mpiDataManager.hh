// Include Guard
#ifndef mpiDataManager_hh
#define mpiDataManager_hh

#ifdef PLB_MPI_PARALLEL
#include <modules/mpiDataManager.h>
#include <parallelism/mpiManager.h>
#include <atomicBlock/dataField3D.hh>
#include <exception>

namespace plb{
	namespace global{

	MpiDataManager::MpiDataManager(){}

	MpiDataManager::~MpiDataManager(){}

	template<>
	void MpiDataManager::sendScalarField3D<int>(const ScalarField3D<int>& field, const Box3D& fromDomain){
		PLB_PRECONDITION( contained(fromDomain, field.getBoundingBox()) );
		const pluint nDataPacks = fromDomain.nCells();
		if(nDataPacks==0){return;}
		// Create mpi variables
		const int count=4;
		const int rank = mpi().getRank();
		plint Xmin = 0; plint Xmax = 0; plint Ymin = 0; plint Ymax = 0; plint Zmin = 0; plint Zmax = 0;
		if(fromDomain.x0 == fromDomain.x1 || fromDomain.y0 == fromDomain.y1 || fromDomain.z0 == fromDomain.z1){
			throw std::domain_error("Domain Boundary Mismatch MpiDataManager::sendScalarFiel3D"); }
		if(fromDomain.x0 < fromDomain.x1){ Xmin = fromDomain.x0; Xmax = fromDomain.x1;}
		else{ Xmin = fromDomain.x1; Xmax = fromDomain.x0; }
		if(fromDomain.y0 < fromDomain.y1){ Ymin = fromDomain.y0; Ymax = fromDomain.y1;}
		else{ Ymin = fromDomain.y1; Ymax = fromDomain.y0; }
		if(fromDomain.z0 < fromDomain.z1){ Zmin = fromDomain.z0; Zmax = fromDomain.z1;}
		else{ Zmin = fromDomain.z1; Zmax = fromDomain.z0; }

		for (plint iX=Xmin; iX<=Xmax; ++iX) {
			for (plint iY=Ymin; iY<=Ymax; ++iY) {
				for (plint iZ=Zmin; iZ<=Zmax; ++iZ) {
					// Create a buffer for the data
					long* sendBuffer = new long[count];
					// Sync mpi processes
					mpi().barrier();
					// Fill the buffer
					sendBuffer[0] = iX; sendBuffer[1]= iY; sendBuffer[2]= iZ; sendBuffer[3]=field.get(iX,iY,iZ);
					// Send the data to all other processes
					mpi().bCast(sendBuffer, count, rank);
					// Clear the buffer
					delete[] sendBuffer;
				}
			}
		}
	}

	template<>
	void MpiDataManager::receiveScalarField3D<int>(ScalarField3D<int>& field, const Box3D& fromDomain, const int& fromId) const{
		PLB_PRECONDITION( contained(fromDomain, field.getBoundingBox()) );
		const int count = 4;
		const int rank = mpi().getRank();
		const pluint nDataPacks = fromDomain.nCells();
		int n = 0;
		while(n < nDataPacks){
			// Create a buffer for the data
			long* recvBuffer = new long[count];
			// Sync the mpi processes
			mpi().barrier();
			// Receive the Data
			mpi().receive(recvBuffer, count, fromId);
			// Put the data into the ScalarField
			if(sizeof(recvBuffer[3]) <= sizeof(int)){
				field.get(recvBuffer[0], recvBuffer[1], recvBuffer[2]) = (int)recvBuffer[3];
			}
			else{ throw std::overflow_error("MpiDataManager::receiveScalarField3D<int>");}
			// clear the Buffer
			delete[] recvBuffer;
			n++;
		}

	}

	} // namespace global
} // namespace plb
#endif
#endif
