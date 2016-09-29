#ifndef CONSTANTS_HH
#define CONSTANTS_HH

#include "constants.h"
#include <palabos3D.hh>
#include "myheaders3D.hh"

namespace plb{

	template<typename T>
	Constants<T>::Constants()
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing Constants";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->parameterXmlFileName = ""; this->epsilon=0; this->minRe = 0; this->maxRe=0;this->maxGridLevel=0;
			this->margin=0; this->borderWidth=0; this->extraLayer=0; this->blockSize=0; this->envelopeWidth=0; this->initialTemperature=0;
			this->gravitationalAcceleration=0; this->master = false;
			this->test = true; this->testRe = 0; this->testTime = 20; this->maxT = 0; this->imageSave = 0; this->testIter = 0;
			this->master = global::mpi().isMainProcessor();
			this->c.reset(this);
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Constants already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T>
	Constants<T>::~Constants()
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Destroying Constants";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
		objCount--;
	}

	template<typename T>
	void Constants<T>::initialize(const std::string& fileName)
	{
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Creating Constants";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			this->parameterXmlFileName = fileName;
			XMLreader r(fileName);
			std::string dir;
			r["dir"]["out"].read(dir);
			global::directories().setOutputDir(dir);
			r["dir"]["log"].read(dir);
			global::directories().setLogOutDir(dir);
			r["dir"]["image"].read(dir);
			global::directories().setImageOutDir(dir);

			r["lbm"]["maxRe"].read(this->maxRe);
			r["lbm"]["minRe"].read(this->minRe);
			r["lbm"]["epsilon"].read(this->epsilon);
			r["lbm"]["gravity"].read(this->gravitationalAcceleration);
			r["lbm"]["u0lb"].read(lb.u);
			r["lbm"]["u0"].read(physical.u);
			r["lbm"]["pl"].read(physical.resolution);

			r["refinement"]["margin"].read(this->margin);
			r["refinement"]["borderWidth"].read(this->borderWidth);
			if(this->margin < this->borderWidth){throw marginEx;} // Requirement: margin>=borderWidth.
			r["refinement"]["maxGridLevel"].read(this->maxGridLevel);
			r["refinement"]["extraLayer"].read(this->extraLayer);
			r["refinement"]["blockSize"].read(this->blockSize);
			r["refinement"]["envelopeWidth"].read(this->envelopeWidth);

			std::string val;
			r["wall"]["meshFileName"].read(val);
			wall.fileName = val;
			r["wall"]["dynamicMesh"].read(wall.dynamicMesh);
			r["wall"]["material"].read(wall.material);
			r["wall"]["referenceDirection"].read(wall.referenceDirection);
			r["wall"]["initialTemperature"].read(wall.initialTemperature);
			Array<T,3> dim = Array<T,3>(0,0,0);
			r["wall"]["lx"].read(dim[0]);
			r["wall"]["ly"].read(dim[1]);
			r["wall"]["lz"].read(dim[2]);
			wall.dim = dim;

			r["obstacle"]["meshFileName"].read(val);
			obstacle.fileName = val;
			r["obstacle"]["dynamicMesh"].read(obstacle.dynamicMesh);
			r["obstacle"]["material"].read(obstacle.material);
			Array<T,3> start = Array<T,3>(0,0,0);
			r["obstacle"]["obstacleStartX"].read(start[0]);
			r["obstacle"]["obstacleStartY"].read(start[1]);
			r["obstacle"]["obstacleStartZ"].read(start[2]);
			obstacle.start = start;
			r["obstacle"]["referenceDirection"].read(obstacle.referenceDirection);
			r["obstacle"]["density"].read(obstacle.density);
			dim = Array<T,3>(0,0,0);
			r["obstacle"]["lx"].read(dim[0]);
			r["obstacle"]["ly"].read(dim[1]);
			r["obstacle"]["lz"].read(dim[2]);
			obstacle.dim = dim;
			T l = 1000;
			for(int i = 0; i<3; i++){
				if(wall.dim[i] < l){ l = wall.dim[i]; }
				if(obstacle.dim[i] < l){ l = obstacle.dim[i]; }
			}
			physical.length = l;
			lb.dx = physical.length / physical.resolution;
			lb.lx = wall.dim[0] / lb.dx;
			physical.lx = wall.dim[0];
			lb.ly = wall.dim[1] / lb.dx;
			physical.ly = wall.dim[1];
			lb.lz = wall.dim[2] / lb.dx;
			physical.lz = wall.dim[2];
			//wall.dim = wall.dim / dx;
			//obstacle.dim = obstacle.dim / dx;

			int tmp = 0;
			r["simulation"]["test"].read(tmp);
			if(tmp == 1){ this->test = true; }else{ this->test = false;}
			r["simulation"]["testRe"].read(this->testRe);
			r["simulation"]["testTime"].read(this->testTime);
			r["simulation"]["maxT"].read(this->maxT);
			r["simulation"]["imageSave"].read(this->imageSave);
			r["simulation"]["testIter"].read(this->testIter);
			r["simulation"]["ibIter"].read(this->ibIter);
			int prec = 0;
			r["simulation"]["precision"].read(prec);
			r["simulation"]["initialTemperature"].read(this->initialTemperature);
			switch(prec){
				case 1: this->precision = Precision::FLT;
				case 2: this->precision = Precision::DBL;
				case 3: this->precision = Precision::LDBL;
			}
			double ratio = 3;
			double x = 1;
			double y = 3;
			double cs = sqrt(x / y);
			double maxMach = 0.1;
			// Fill the 2D array with standard values
			for(plint grid = 0; grid <= this->maxGridLevel; grid++){
				try{
					T resolution = physical.resolution * util::twoToThePowerPlint(grid);
					T scaled_u0lb = lb.u / util::twoToThePowerPlint(grid);
					double mach = scaled_u0lb / cs;
					if(mach > maxMach){std::cout<<"Local Mach= "<<mach<<"\n"; throw localMachEx;}
					if(resolution == 0){throw resolEx;}
					if(this->test){
						this->minRe = this->testRe; this->maxRe = this->testRe+1;
						IncomprFlowParam<T> p = IncomprFlowParam<T>(physical.u,scaled_u0lb,testRe,physical.length,resolution,lb.lx,lb.ly,lb.lz);
						// Check local speed of sound constraint
						T dt = p.getDeltaT();
						T dx = p.getDeltaX();
						if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
					}
					else{
						for(int reynolds = 0; reynolds <= this->maxRe; reynolds++){
							IncomprFlowParam<T> p =
								IncomprFlowParam<T>(physical.u,scaled_u0lb,reynolds,physical.length,resolution,lb.lx,lb.ly,lb.lz);
							// Check local speed of sound constraint
							T dt = p.getDeltaT();
							T dx = p.getDeltaX();
							if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
						}
					}
				}
				catch(const std::exception& e){
					if(grid == 0){ throw e; }
					else{ this->maxGridLevel = grid - 1; break; }
				}
			}
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] Done Creating Constants";
				if(this->master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			objCount++;
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
	}


} // namespace plb

#endif // CONSTANTS_HH
