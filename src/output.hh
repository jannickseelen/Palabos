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
#include <iostream>
#include <iomanip>

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
		try{
			bool parallel = false;
			#ifdef PLB_MPI_PARALLEL
				parallel = true;
			#endif
			global::log().init(parallel);
			master = global::mpi().isMainProcessor();
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::elapsedTime()
	{
		try{
			if(Constants<T>::test){std::thread(timeLoop);}
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::timeLoop()
	{
		try{
			double elapsed = 0;
			double testTime = Constants<T>::testTime*60;
			bool running = true;
			while(running){
				elapsed = timer.getTime() - startTime;
				if(elapsed > testTime){ throw std::runtime_error("Test Time Limit Exceeded!"); running = false; }
				std::this_thread::sleep_for(std::chrono::minutes(1));
			}
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::writeGif()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Writing GIF";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			const plint imSize = 600;
			plint nx = Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getNx();
			plint ny = Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getNy();
			plint nz = Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice->getNz();
			nx = nx/2;
			ny = ny/2;
			nz = nz/2;
			Box3D slice(nx, nx-1, ny, ny-1, nz, nz);
			ImageWriter<T> imageWriter("leeloo");
			imageWriter.writeScaledGif(createFileName("u", gifCount, 6),
				*computeVelocityNorm(*Variables<T,BoundaryType,SurfaceData,Descriptor>::lattice, slice),imSize, imSize );
			gifCount++;
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Writing GIF";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}
	
	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::writeImages()
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Writing VTK";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			T dx = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaX();
			T dt = Variables<T,BoundaryType,SurfaceData,Descriptor>::p.getDeltaT();
			// Create a VTK data object and indicate the cell width through the
			//   parameter dx. This object is of the same type as the simulation, type T.
			VtkImageOutput3D<T,float,3> vtkOut("simulationData.dat", dx);

			// Add a 3D tensor-field to the VTK file, representing the velocity, and rescale
			//   with the units of a velocity, dx/dt. Explicitly convert the data to single-
			//   precision floats in order to save storage space.
			int size = Variables<T,BoundaryType,SurfaceData,Descriptor>::velocity.size();
			plb_ofstream v_file("velocity.dat");
			for( int i = 0; i<size; i++){
				vtkOut.writeData(Variables<T,BoundaryType,SurfaceData,Descriptor>::velocity[i], "velocity", dx/dt);
				v_file << std::setprecision(10) << Variables<T,BoundaryType,SurfaceData,Descriptor>::velocity[i] << std::endl;
			}

			// Add another 3D tensor-field for the vorticity, again as floats.
			size = Variables<T,BoundaryType,SurfaceData,Descriptor>::vorticity.size();
			plb_ofstream w_file("vorticity.dat");
			for( int i = 0; i<size; i++){
				vtkOut.writeData(Variables<T,BoundaryType,SurfaceData,Descriptor>::vorticity[i],	"vorticity", 1./dt);
				w_file << std::setprecision(10) << Variables<T,BoundaryType,SurfaceData,Descriptor>::vorticity[i] << std::endl;
			}

			// To end-with add a scalar-field for the density.
			size = Variables<T,BoundaryType,SurfaceData,Descriptor>::density.size();
			plb_ofstream r_file("density.dat");
			for( int i = 0; i<size; i++){
				vtkOut.writeData(Variables<T,BoundaryType,SurfaceData,Descriptor>::density[i], "density", 1.);
				r_file << std::setprecision(10) << Variables<T,BoundaryType,SurfaceData,Descriptor>::density[i] << std::endl;
			}
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Writing VTK";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::startMessage()
	{
		try
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
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::simMessage()
	{
		try{
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
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

	template<typename T, class BoundaryType, class SurfaceData, template<class U> class Descriptor>
	void Output<T,BoundaryType,SurfaceData,Descriptor>::stopMessage()
	{
		try{
		std::string mesg = "SIMULATION COMPLETE";
		if(master){std::cout << mesg << std::endl;}
		global::log(mesg);
		elapsedTime();
		timer.stop();
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}

} // namespace plb

#endif // OUTPUT_HH

