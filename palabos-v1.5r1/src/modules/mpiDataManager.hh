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

	void MpiDataManager::checkDomain(int rank, Box3D domain){
		if((domain.x1 <= domain.x0) || (domain.y1 <= domain.y0) || (domain.z1 <= domain.z0)){
			std::string err_str("Rank= ");
			err_str.append(std::to_string(rank));
			err_str.append(" Domain=[");
			err_str.append(std::to_string(domain.x0));
			err_str.append(",");
			err_str.append(std::to_string(domain.x1));
			err_str.append("][");
			err_str.append(std::to_string(domain.y0));
			err_str.append(",");
			err_str.append(std::to_string(domain.y1));
			err_str.append("][");
			err_str.append(std::to_string(domain.z0));
			err_str.append(",");
			err_str.append(std::to_string(domain.z1));
			err_str.append("]");
			throw std::runtime_error(err_str);
		}
	}

	template<>
	void MpiDataManager::sendScalarField3D<int>(const ScalarField3D<int>& field, const Box3D& fromDomain){
		try{
			PLB_PRECONDITION( contained(fromDomain, field.getBoundingBox()) );
			const pluint nDataPacks = fromDomain.nCells();
			if(nDataPacks==0){return;}
			// Create mpi variables
			const int count=4;
			const int rank = mpi().getRank();
			#ifdef PLB_DEBUG
				std::cout << "Rank " << rank << " Sending " << nDataPacks << " DataPacks " << std::endl;
			#endif
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
						// mpi().barrier();
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
		catch(const std::exception& e){
			std::cout << "Excpetion Caught in MpiDataManager::sendScalarField3D<int>:" << e.what() << std::endl;
			throw e;
		}
	}

	template<>
	void MpiDataManager::receiveScalarField3D<int>(ScalarField3D<int>& field, const Box3D& fromDomain, const int& fromId) const{
		PLB_PRECONDITION( contained(fromDomain, field.getBoundingBox()) );
		const int count = 4;
		const int rank = mpi().getRank();
		const pluint nDataPacks = fromDomain.nCells();
		int n = 0;
		#ifdef PLB_DEBUG
			std::cout << "Rank " << rank << " Receiving " << nDataPacks << "DataPacks from Rank " << fromId << std::endl;
		#endif
		while(n < nDataPacks){
			// Create a buffer for the data
			long* recvBuffer = new long[count];
			// Sync the mpi processes
			// mpi().barrier();
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

	std::vector<Box3D> MpiDataManager::splitDomains(const Box3D& domain){
		int minX, maxX, minY, maxY, minZ, maxZ;
		if(domain.x0 < 0){ minX = 0; }else{if(domain.x0 < domain.x1){ minX = domain.x0;}else{minX = domain.x1;}}
		if(domain.x1 >= domain.getNx()){maxX = domain.getNx();}else{if(domain.x1 > domain.x0){maxX = domain.x1;}else{maxX = domain.x0;}}
		if(domain.y0 < 0){ minY = 0; }else{if(domain.y0 < domain.y1){minY = domain.y0;}else{ minY = domain.y1;}}
		if(domain.y1 >= domain.getNy()){maxY = domain.getNy();}else{if(domain.y1 > domain.y0){maxY = domain.y1;}else{maxY = domain.y0;}}
		if(domain.z0 < 0){ minZ = 0; }else{if(domain.z0 < domain.z1){minZ = domain.z0;}else{minZ = domain.z1;}}
		if(domain.z1 >= domain.getNz()){maxZ = domain.getNz();}else{if(domain.z1 > domain.z0){maxZ = domain.z1;}else{maxZ = domain.z0;}}
		const int nproc = mpi().getSize();
		std::vector<Box3D> mpiDomains;
		// mpiDomains.resize(nproc);
		int nSide = std::cbrt(nproc);
		if(nSide == 0){ throw std::runtime_error("Qubic Root of nprocs failed");}
		plint xdif, xrem, ydif, yrem, zdif, zrem;
		bool mpiExcess = false;
		if(maxX - minX < nSide){ xdif = 1; mpiExcess = true; }
		else{ xdif = floor((maxX-minX)/nSide); }
		if(maxY - minY < nSide){ ydif = 1; mpiExcess = true; }
		else{ ydif = floor((maxY-minY)/nSide); }
		if(maxZ - minZ < nSide){ zdif = 1; mpiExcess = true; }
		else{ zdif = floor((maxZ-minZ)/nSide); }
		if((maxX-minX) % nSide !=0){ xrem = (maxX-minX) % nSide; }else{xrem = 0;}
		if((maxY-minY) % nSide !=0){ yrem = (maxY-minY) % nSide; }else{yrem = 0;}
		if((maxZ-minZ) % nSide !=0){ zrem = (maxZ-minZ) % nSide; }else{zrem = 0;}
		if(xdif==0 || ydif==0 || zdif==0){
			std::cout << "Rank " << global::mpi().getRank() << " Domain [" << minX << "," << maxX << "," << minY << "," << maxY << ","
			<< minZ << "," << maxZ << "]" << std::endl;
			throw std::runtime_error("Domain Boundaries are not set properly");
		}
		bool error = false;
		std::vector<std::string> error_domains;
		int r = 0; int y_last = 0; int x_last = 0; int z_last = 0;
		for(int x = 0; x<nSide; x++){
			int x0, x1;
			if(x==0){x0 = minX; x1 = minX + xdif; }
			else{ x0 = x_last + 1; if(x==nSide-1){x1 = maxX;} else{ x1 = x0 + xdif; } }
			for(int y=0; y<nSide; y++){
				int y0, y1;
				if(y==0){y0 = minY; y1 = minY + ydif;}
				else{ y0 = y_last + 1; if(y==nSide-1){y1 = maxY;} else{ y1 = y0 + ydif;} }
				for(int z=0; z<nSide; z++){
					int z0, z1;
					if(z==0){z0 = minZ; z1 = minZ + zdif; }
					else{ z0 = z_last + 1; if(z==nSide-1){z1 = maxZ;} else{ z1 = z0 + zdif; } }
					x_last = x1; y_last = y1; z_last = z1;
					Box3D rankDomain(x0,x1,y0,y1,z0,z1);
					checkDomain(r, rankDomain);
					mpiDomains.push_back(rankDomain);
					r++;
				}
			}
		}
		if(mpiExcess){
			int missing = nproc - mpiDomains.size();
			for(int i = 0; i < missing; i++){
				Box3D empty(0,0,0,0,0,0);
				mpiDomains.push_back(empty);
			}
		}
		if(mpiDomains.size() != nproc){ throw std::runtime_error("One or more mpi Domains where not initialized."); }
		return mpiDomains;
	}

	} // namespace global
} // namespace plb
#endif
#endif
