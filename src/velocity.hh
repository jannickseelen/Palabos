#ifndef VELOCITY_HH
#define VELOCITY_HH

namespace plb{

	template<typename T>
	SurfaceVelocity<T>::SurfaceVelocity(const Array<T,3>& start, const T& mass_, const T& g_)
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			mass = mass_;
			g = -g_;
			location.push_back(start);
			acceleration.push_back(Array<T,3>(0, 0, g));
			force.push_back(Array<T,3>(0,0,g*mass));
			velocity.push_back(Array<T,3>(0,0,0));
			time.push_back((T)0);
			objCount++;
		}
		else
		{
			std::string ex = "Static Class Surface Velocity already defined";
			std::string line = std::to_string(__LINE__);
			std::string mesg = "[ERROR]: "+ex+" [FILE:"+__FILE__+",LINE:"+line+"]";
			global::log(mesg);
			throw std::runtime_error(mesg);
		}
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::update(const T& timeLB, Array<T,3> force){
		Array<T,3> ds;
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			force[2] += g*mass;
			acceleration = force/mass;
			T dt = timeLB - time.back();
			time.push_back(timeLB);
			Array<T,3> v = velocity.back() + acceleration*dt;
			velocity.push_back(v);
			Array<T,3> c = location.back() + velocity*dt;
			ds = c - location.back();
			location.push_back(c);
			#ifdef PLB_DEBUG
				mesg = "[DEBUG] DONE SurfaceVelocity= ["+std::to_string(velocity[0])+","+std::to_string(velocity[1])+
				","+std::to_string(velocity[2])+"]";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return ds;
	}



} // NAMESPACE PLB

#endif // VELOCITY_HH
