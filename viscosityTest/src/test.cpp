//
//  test.cpp
//  Palabos
//
//  Created by MAC on 10/02/16.
//  Copyright © 2016 TUDelft. All rights reserved.
//
//
#ifdef PLB_DEBUG
	#define NDEBUG
#endif

typedef double T;

// PALABOS INCLUDES
#include <latticeBoltzmann/nearestNeighborLattices3D.hh>

// MY INCLUDES
#include "myheaders3D.hh"

//C++ includes
#include <string>

#define Descriptor plb::descriptors::D3Q19Descriptor

typedef plb::Array<T,3> BoundaryType;
typedef plb::Array<T,3> SurfaceData;

int main(int argc, char* argv[])
{
	try{
		//plb::installSigHandler();
		plb::plbInit(&argc, &argv, true); // Initialize Palabos
		bool master = plb::global::mpi().isMainProcessor();
		std::string fileName = "";
		plb::global::argv(argc-1).read(fileName);
		plb::Helper<T,BoundaryType,SurfaceData,Descriptor> h;
		h.initialize(fileName);
		plb::Constants<T>* constants = plb::Constants<T>::c.get();
		plb::Wall<T,BoundaryType,SurfaceData,Descriptor>* wall = plb::Wall<T,BoundaryType,SurfaceData,Descriptor>::w.get();
		plb::Obstacle<T,BoundaryType,SurfaceData,Descriptor>* obstacle = plb::Obstacle<T,BoundaryType,SurfaceData,Descriptor>::o.get();
		plb::Variables<T,BoundaryType,SurfaceData,Descriptor>* variables = plb::Variables<T,BoundaryType,SurfaceData,Descriptor>::v.get();
		plb::Output<T,BoundaryType,SurfaceData,Descriptor>* output = plb::Output<T,BoundaryType,SurfaceData,Descriptor>::out.get();
		output->startMessage();
		output->elapsedTime(); // Initialize Test Timer
		plb::pcout << "Min Reynolds = "<<constants->minRe<<" Max Reynolds = "<<constants->maxRe << std::endl;
		plb::pcout << "Min Grid Level = 1 Max Grid Level = "<<constants->maxGridLevel << std::endl;
		for(plb::plint reynolds = constants->minRe; reynolds <= constants->maxRe; reynolds++){
			for(plb::plint gridLevel = 1; gridLevel<= constants->maxGridLevel; gridLevel++){

				variables->update(gridLevel,reynolds);
				variables->setLattice();
				bool converged = false;
				for(int i=0; converged == false; i++)
				{
					variables->iter++;
					variables->time = i + 1.0;
					variables->lattice->collideAndStream();
					variables->updateLattice();
					//if(variables->checkConvergence()){ converged = true; break; }
					#ifdef PLB_DEBUG
						std::string mesg="N collisions="+std::to_string(variables->iter);
						if(master){std::cout << mesg << std::endl;}
						plb::global::log(mesg);
						if(master){std::cout<<"Grid Level="+std::to_string(gridLevel);}
						if(master){std::cout << mesg << std::endl;}
						plb::global::log(mesg);
					#endif
					if(constants->test){ if(variables->iter > constants->testIter){ break; }}
				}
				output->writeImages();
			}
			if(constants->test){ break; }
		}
		output->stopMessage();
		return 0;																	// Return Process Completed
	}
	catch(const std::exception& e){plb::exHandler(e,__FILE__,__FUNCTION__,__LINE__);	return -1;	}
};
