// Include Guard
#ifndef mpiDataManager_hh
#define mpiDataManager_hh

#ifdef PLB_MPI_PARALLEL
#include <modules/mpiDataManager.h>
#include <parallelism/mpiManager.h>
#include <atomicBlock/dataField3D.hh>
#include <string>
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

	std::vector<Box3D> MpiDataManager::splitDomains(const Box3D& domain, const plint& minX, const plint& maxX, const plint& minY, const plint& maxY,
			const plint& minZ, const plint& maxZ){
		const int nproc = mpi().getSize();
		const int rank = mpi().getRank();
		std::vector<Box3D> mpiDomains;
		// mpiDomains.resize(nproc);
		int nSide = std::cbrt(nproc);
		std::vector<int> dx, dy, dz;
		dx.resize(nSide); dy.resize(nSide); dz.resize(nSide);
		plint xdif = (maxX-minX);
		plint ydif = (maxY-minY);
		plint zdif = (maxZ-minZ);
		if(xdif/nSide % 2){
			int x = floor(xdif/nSide);
			for(int n = 0; n!=nSide; n++){if(n==0){dx[n] = x+(xdif/nSide % 2);} else{dx[n] = x;} }
		}
		else{
			int x = floor(xdif/nSide);
			for(int n = 0; n!=nSide; n++){ dx[n] = x; }
		}
		if(ydif/nSide % 2){
			int y = floor(ydif/nSide);
			for(int n = 0; n!=nSide; n++){if(n==0){dy[n] = y+(ydif/nSide % 2);} else{dy[n] = y;} }
		}
		else{
			int y = floor(ydif/nSide);
			for(int n = 0; n!=nSide; n++){ dy[n] = y; }
		}
		if(zdif/nSide % 2){
			int z = floor(zdif/nSide);
			for(int n = 0; n!=nSide; n++){if(n==0){dz[n] = z+(zdif/nSide % 2);} else{dz[n] = z;} }
		}
		else{
			int z = floor(zdif/nSide);
			for(int n = 0; n!=nSide; n++){ dz[n] = z; }
		}
		for(int n = 0; n!=nSide; n++){
			if(dx[n]==0 || dy[n]==0 || dz[n]==0){ throw std::runtime_error("Domain Boundaries are not set properly");}
		}
		bool error = false;
		std::vector<std::string> error_domains;
		int r, y_last, x_last, z_last;
		for(int x = 0; x<nSide; x++){
			int x0, x1;
			if(x==0){x0 = minX;} else{ x0 = x_last + 1; } // Ensures no overlap
			if(x==nSide-1){x1 = maxX;} else{ x1 = x0 + dx[x]; }
			for(int y=0; y<nSide; y++){
				int y0, y1;
				if(y==0){y0 = minY;} else{ y0 = y_last + 1; }
				if(y==nSide-1){y1 = maxY;} else{ y1 = y0 + dy[y];}
				for(int z=0; z<nSide; z++){
					int z0, z1;
					if(z==0){z0 = minZ;} else{ z0 = z_last + 1; }
					if(z==nSide-1){z1 = maxZ;} else{ z1 = z0 + dz[z]; }
					x_last = x1; y_last = y1; z_last = z1;
					Box3D rankDomain(x0,x1,y0,y1,z0,z1);
					if((x1 <= x0) || (y1 <= y0) || (z1 <= z0)){
						error = true;
						std::string err_str("Rank= ");
						err_str.append(std::to_string(r));
						err_str.append(" Domain=[");
						err_str.append(std::to_string(x0));
						err_str.append(",");
						err_str.append(std::to_string(x1));
						err_str.append("][");
						err_str.append(std::to_string(y0));
						err_str.append(",");
						err_str.append(std::to_string(y1));
						err_str.append("][");
						err_str.append(std::to_string(z0));
						err_str.append(",");
						err_str.append(std::to_string(z1));
						err_str.append("]");
						error_domains.push_back(err_str);
					}
					else{mpiDomains.push_back(rankDomain);}
					r++;
				}
			}
		}
		if(error){
			std::string err_str;
			for(int i = 0; i!=error_domains.size(); i++){
				err_str.append(error_domains[i]);
				if(i<error_domains.size()-1){err_str.append(" \n");}
			}
			throw std::runtime_error(err_str);
		}
		if(mpiDomains.size() != nproc){ throw std::runtime_error("One or more mpi Domains where not initialized."); }
		return mpiDomains;
	}

	} // namespace global
} // namespace plb
#endif
#endif
