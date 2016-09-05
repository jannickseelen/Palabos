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
			this->parameterXmlFileName = ""; this->u0lb=0; this->epsilon=0; this->minRe = 0; this->maxRe=0;this->maxGridLevel=0;
			this->referenceResolution=0;
			this->margin=0; this->borderWidth=0; this->extraLayer=0; this->blockSize=0; this->envelopeWidth=0; this->initialTemperature=0;
			this->gravitationalAcceleration=0; this->dynamicWall = false; this->dynamicObstacle = false; this->master = false;
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
			r["lbm"]["u0lb"].read(this->u0lb);
			double x = 1;
			double y = 3;
			double cs = sqrt(x / y);
			double mach = this->u0lb / cs;
			double maxMach = 0.1;
			if(mach > maxMach){throw machEx;}
			r["lbm"]["maxRe"].read(this->maxRe);
			r["lbm"]["minRe"].read(this->minRe);
			r["lbm"]["epsilon"].read(this->epsilon);
			r["simulation"]["initialTemperature"].read(this->initialTemperature);

			r["wall"]["meshFileName"].read(wall_file);
			r["wall"]["dynamicMesh"].read(this->dynamicWall);
			r["wall"]["material"].read(wall_mat);
			r["wall"]["referenceDirection"].read(wall_data[0]);
			r["wall"]["initialTemperature"].read(wall_data[1]);

			r["obstacle"]["meshFileName"].read(obstacle_file);
			r["obstacle"]["dynamicMesh"].read(this->dynamicObstacle);
			r["obstacle"]["material"].read(obstacle_mat);
			r["obstacle"]["obstacleStartX"].read(obstacle_data[0]);
			r["obstacle"]["obstacleStartY"].read(obstacle_data[1]);
			r["obstacle"]["obstacleStartZ"].read(obstacle_data[2]);
			r["obstacle"]["referenceDirection"].read(obstacle_data[3]);
			r["obstacle"]["density"].read(obstacle_data[4]);

			r["refinement"]["margin"].read(this->margin);
			r["refinement"]["borderWidth"].read(this->borderWidth);
			if(this->margin < this->borderWidth){throw marginEx;} // Requirement: margin>=borderWidth.
			r["refinement"]["maxGridLevel"].read(this->maxGridLevel);
			r["refinement"]["extraLayer"].read(this->extraLayer);
			r["refinement"]["blockSize"].read(this->blockSize);
			r["refinement"]["envelopeWidth"].read(this->envelopeWidth);
			r["refinement"]["referenceResolution"].read(this->referenceResolution);
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
			switch(prec){
				case 1: this->precision = Precision::FLT;
				case 2: this->precision = Precision::DBL;
				case 3: this->precision = Precision::LDBL;
			}
			if(this->test){ this->minRe = this->testRe; this->maxRe = this->testRe+1; this->maxGridLevel = 0; }
			double ratio = 3;
			// Fill the 2D array with standard values
			for(plint grid = 0; grid <= this->maxGridLevel; grid++){
				T resolution = this->referenceResolution * util::twoToThePowerPlint(grid);
				T scaled_u0lb = this->u0lb * util::twoToThePowerPlint(grid);
				mach = scaled_u0lb / cs;
				if(mach > maxMach){std::cout<<"Local Mach= "<<mach<<"\n"; throw localMachEx;}
				if(resolution == 0){throw resolEx;}
				if(this->test){
					IncomprFlowParam<T> p = IncomprFlowParam<T>(scaled_u0lb,testRe,resolution,1,1,1);
					// Check local speed of sound constraint
					T dt = p.getDeltaT();
					T dx = p.getDeltaX();
					if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
				}
				else{
					for(int reynolds = 0; reynolds <= this->maxRe; reynolds++){
						IncomprFlowParam<T> p = IncomprFlowParam<T>(scaled_u0lb,reynolds,resolution,1,1,1);
						// Check local speed of sound constraint
						T dt = p.getDeltaT();
						T dx = p.getDeltaX();
						if(dt > (dx / sqrt(ratio))){std::cout<<"dt:"<<dt<<"<(dx:"<<dx<<"/sqrt("<<ratio<<")"<<"\n"; throw superEx;}
					}
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
