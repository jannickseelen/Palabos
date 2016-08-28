#ifndef OUTPUT_HH
#define OUTPUT_HH

#include "output.h"
#include <palabos3D.hh>
#include "myheaders3D.hh"

#include <thread>
#include <ctime>
#include <chrono>
#include <string>
#include <exception>

namespace plb{

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Output<T,BoundaryType,SurfaceData,Descriptor>::Output()
	{
		if(objCount==0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing Output";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->out.reset(this);
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Output already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	Output<T,BoundaryType,SurfaceData,Descriptor>::~Output()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Output";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		objCount--;
	}


	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::initialize()
	{
		bool parallel = false;
		#ifdef PLB_MPI_PARALLEL
			parallel = true;
		#endif
		global::log().init(parallel);
		master = global::mpi().isMainProcessor();
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::elapsedTime()
	{
		if(Constants<T>::test){std::thread(timeLoop);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::timeLoop()
	{
		double elapsed = 0;
		double testTime = Constants<T>::testTime*60;
		bool running = true;
		while(running){
			elapsed = timer.getTime() - startTime;
			if(elapsed > testTime){ throw std::runtime_error("Test Time Limit Exceeded!"); running = false; }
			std::this_thread::sleep_for(std::chrono::minutes(1));
		}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::writeGifs(MultiBlockLattice3D<T,Descriptor>& lattice, plint iter)
	{
		const plint imSize = 600;
		const plint nx = lattice.getNx();
		const plint ny = lattice.getNy();
		const plint nz = lattice.getNz();
		Box3D slice(0, nx-1, 0, ny-1, nz/2, nz/2);
		ImageWriter<T> imageWriter("leeloo");
		imageWriter.writeScaledGif(createFileName("u", iter, 6),*computeVelocityNorm(lattice, slice),imSize, imSize );
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::writeImages()
	{
		/*T dx = Constants<T>::parameters[v->gridLevel][v->reynolds]->getDeltaX();
		T dt = Constants<T>::parameters[v->gridLevel][v->reynolds]->getDeltaT();
		int n = v->velocity.size();
		for(int i=0; i<n; i++){
			VtkImageOutput3D<T> vtkOut(createFileName("vtk",i,n), dx);
			MultiTensorField3D<T,3> v_field = v->velocity[i];
			vtkOut.writeData<3,float>(v_field, name, dx/dt);
		}*/
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::startMessage()
	{
		time_t rawtime;
		struct tm * timeinfo;
		time (&rawtime);
		timeinfo = localtime (&rawtime);
		if(master)
		{
			std::cout<<"SIMULATION START "<< asctime(timeinfo) << std::endl;
		}
		#ifdef PLB_MPI_PARALLEL
		if(master){
			int size = plb::global::mpi().getSize();
			std::cout<<"NUMBER OF MPI PROCESSORS="<< size << std::endl;
			//if((int)cbrt(size) % 1){ throw std::runtime_error("Number of MPI Processess must satisfy Cubic Root");}
			std::string imaster =  master ? " YES " : " NO ";
			std::cout<<"Is this the main process?"<< imaster << std::endl;
		}
		#endif
		std::string outputDir = "./tmp/";
		global::directories().setOutputDir(outputDir);// Set output DIR w.r.t. to current DIR
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::simMessage()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Creating Timer";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		timer.start();		// Start Timer
		startTime = timer.getTime();
		#ifdef PLB_DEBUG
			mesg = "[DEBUG] Timer Started";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
			if(Constants<T>::test){mesg = "[DEBUG] Starting Test";}
			else{mesg = "[DEBUG] Starting Normal Run";}
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::stopMessage()
	{
		std::string mesg = "SIMULATION COMPLETE";
		if(master){std::cout << mesg << std::endl;}
		global::log(mesg);
		elapsedTime();
		timer.stop();
	}

} // namespace plb

#endif // OUTPUT_HH

