#ifndef VELOCITY_HH
#define VELOCITY_HH

namespace plb{

	template<typename T>
	SurfaceVelocity<T>::SurfaceVelocity()
	{
		if(objCount == 0)
		{
			master = global::mpi().isMainProcessor();
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Constructing SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
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
	void SurfaceVelocity<T>::initialize(const Array<T,3>& start, const T& mass_, const T& g_)
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Initializing SurfaceVelocity";
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
		#ifdef PLB_DEBUG
			mesg = "[DEBUG] DONE Initializing SurfaceVelocity";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::update(const T& timeLB, Array<T,3> force){
		Array<T,3> ds = Array<T,3>();
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			force[2] += g*mass;
			T dt = timeLB - time.back();
			time.push_back(timeLB);
			Array<T,3> a = Array<T,3>();
			Array<T,3> v = Array<T,3>();
			Array<T,3> c = Array<T,3>();
			for(int i = 0; i<3; i++){
				a[i] = force[i] / mass;
				v[i] = velocity.back()[i] + a[i]*dt;
				c[i] = location.back()[i] + v[i]*dt;
				ds[i] = c[i] - location.back()[i];
			}
			acceleration.push_back(a);
			velocity.push_back(v);
			location.push_back(c);
			#ifdef PLB_DEBUG
				T v0 = v[0];
				T v1 = v[1];
				T v2 = v[2];
				mesg = "[DEBUG] DONE SurfaceVelocity= ["+std::to_string(v0)+","+std::to_string(v1)+
				","+std::to_string(v2)+"]";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return ds;
	}



} // NAMESPACE PLB

#endif // VELOCITY_HH
