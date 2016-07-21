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
			const int nprocs = mpi().getSize();
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
						for(int n = 0; n<nprocs; n++){ if(n != rank){ mpi().send(sendBuffer,4,n); } }
						// Clear the buffer
						delete[] sendBuffer;
					}
				}
			}
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}

	template<>
	void MpiDataManager::receiveScalarField3D<int>(ScalarField3D<int>& field, const Box3D& fromDomain, const int& fromId) const{
		try{
			PLB_PRECONDITION( contained(fromDomain, field.getBoundingBox()) );
			const int count = 4;
			const int rank = mpi().getRank();
			const pluint nDataPacks = fromDomain.nCells();
			int n = 0;
			#ifdef PLB_DEBUG
				std::cout << "Rank " << rank << " Receiving " << nDataPacks << " DataPacks from Rank " << fromId << std::endl;
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
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}

	template<>
	TriangleSet<double> MpiDataManager::receiveTriangleSet<double>(){
		TriangleSet<double> triangles;
		try{
			if(!mpi().isMainProcessor()){
				const int rank = mpi().getRank();
				if(rank == 0){ throw std::runtime_error("Error in receiveTriangleSet, Master should not receive");}
				#ifdef PLB_DEBUG
					std::cout << "[DEBUG] Rank "<<rank<<" is trying to receive a triangleSet"<<std::endl;
				#endif
				std::vector<Array<Array<double,3>,3> > list;
				plint* buffer = new plint[1];
				mpi().receive(buffer, 1, 0);
				plint nTriangles = *buffer;
				delete[] buffer;
				for(int i=0; i<nTriangles; i++){
					Array<Array<double,3>,3> triangle;
					for(int p =0; p<3; p++){
						double* recvBuffer = new double[3];
						Array<double, 3> point;
						mpi().receive(recvBuffer,3,0);
						point[0] = recvBuffer[0]; point[1] = recvBuffer[1]; point[2] = recvBuffer[2];
						triangle[p] = point;
						delete[] recvBuffer;
					}
					list.push_back(triangle);
				}
				triangles = TriangleSet<double>(list, DBL);
				#ifdef PLB_DEBUG
					std::cout << "[DEBUG] Rank "<<rank<<" received a triangleSet with "<<list.size()<<" Triangles"<<std::endl;
				#endif
			}
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
		return triangles;
	}

	template<>
	void MpiDataManager::sendTriangleSet(const TriangleSet<double>& triangles){
		try{
			if(mpi().isMainProcessor()){
				const int rank = mpi().getRank();
				if(rank != 0){ throw std::runtime_error("Error in sendTriangleSet, process is not Master");}
				std::vector<Array<Array<double,3>,3> > list = triangles.getTriangles();
				plint nTriangles = list.size();
				#ifdef PLB_DEBUG
					std::cout << "[DEBUG] Master is sending a triangleSet with "<< nTriangles<<" Triangles"<< std::endl;
				#endif
				const int nprocs = mpi().getSize();
				plint* buffer = new plint[1];
				buffer[0] = nTriangles;
				for(int n = 0; n<nprocs; n++){ if(n != rank){ mpi().send(buffer,1,n); } }
				delete[] buffer;
				for(int i=0; i<nTriangles; i++){
					Array<Array<double,3>,3> triangle = list[i];
					for(int p = 0; p<3; p++){
						double* sendBuffer = new double[3];
						Array<double,3> point = triangle[p];
						sendBuffer[0] = point[0]; sendBuffer[1]=point[1]; sendBuffer[2]=point[2];
						for(int n = 0; n<nprocs; n++){ if(n != rank){ mpi().send(sendBuffer,3,n); } }
						delete[] sendBuffer;
					}
				}
			}
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
	}

	Array<int,2> MpiDataManager::checkIfCrossed(const int& fMin, const int& fMax, const int& n){
		Array<int,2> boundary;
		try{
			if(fMin == fMax){ throw std::runtime_error("Exception in MpiDataManager::checkIfCrossed domain limits are equal");}
			int min = 0; int max=0;
			if(fMin<0){ min = 0; }
			else{
				if(fMin<fMax){
					min = fMin;
				}
				else{
					min = fMax;
				}
			}
			if(fMax >= n){
				max= n;
			}
			else{
				if(fMax > fMin){
					max = fMax;
				}
				else{
					max = fMin;
				}
			}
			boundary[0] = min;
			boundary[1] = max;
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
		return boundary;
	}

	std::vector<Box3D> MpiDataManager::splitDomains(const Box3D& domain){
		std::vector<Box3D> mpiDomains;
		try{
			int minX, maxX, minY, maxY, minZ, maxZ;
			Array<int,2> boundary;
			boundary = checkIfCrossed(domain.x0, domain.x1, domain.getNx());
			minX = boundary[0]; maxX = boundary[1];
			boundary = checkIfCrossed(domain.y0, domain.y1, domain.getNy());
			minY = boundary[0]; maxY = boundary[1];
			boundary = checkIfCrossed(domain.z0, domain.z1, domain.getNz());
			minZ = boundary[0]; maxZ = boundary[1];
			const int nproc = mpi().getSize();
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
			if(!mpiExcess){ mpiDomains.resize(nproc);}
			if(xdif==0 || ydif==0 || zdif==0){
				std::cout << "Rank " << global::mpi().getRank() << " Domain [" << minX << "," << maxX << "," << minY << "," << maxY << ","
				<< minZ << "," << maxZ << "]" << std::endl;
				throw std::runtime_error("Domain Boundaries are not set properly");
			}
			bool error = false;
			std::vector<std::string> error_domains;
			int r = 0; int y_last = 0; int x_last = 0; int z_last = 0;
			for(int x = 0; x<nSide; x++){
				if(r == nproc){ break;}
				if(x_last == maxX){ break; }
				int x0, x1;
				if(x==0){x0 = minX; x1 = minX + xdif; }
				else{ x0 = x_last + 1; if(x==nSide-1){x1 = maxX;} else{ x1 = x0 + xdif; } }
				for(int y=0; y<nSide; y++){
					if(y_last == maxY){ break; }
					int y0, y1;
					if(y==0){y0 = minY; y1 = minY + ydif;}
					else{ y0 = y_last + 1; if(y==nSide-1){y1 = maxY;} else{ y1 = y0 + ydif;} }
					for(int z=0; z<nSide; z++){
						if(z_last == maxZ){ break; }
						int z0, z1;
						if(z==0){z0 = minZ; z1 = minZ + zdif; }
						else{ z0 = z_last + 1; if(z==nSide-1){z1 = maxZ;} else{ z1 = z0 + zdif; } }
						x_last = x1; y_last = y1; z_last = z1;
						if(r == nproc){ break;}
						Box3D rankDomain(x0,x1,y0,y1,z0,z1);
						checkDomain(r, rankDomain);
						if(mpiExcess){mpiDomains.push_back(rankDomain);}
						else{mpiDomains[r] = rankDomain;}
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
		}
		catch(const std::exception& e){
			std::string ex = e.what();
			std::string line = std::to_string(__LINE__);
			global::log().entry("[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]");
			throw e;
		}
		return mpiDomains;
	}

	} // namespace global
} // namespace plb
#endif
#endif
