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
		Array<T,3> a = Array<T,3>(0, 0, g);
		acceleration.push_back(a);
		Array<T,3> f = Array<T,3>(0,0,g*mass);
		force.push_back(f);
		Array<T,3> v = Array<T,3>(0,0,0);
		velocity.push_back(v);
		T t = 0;
		time.push_back(t);
		#ifdef PLB_DEBUG
			mesg = "[DEBUG] DONE Initializing SurfaceVelocity";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::update(const T& timeLB, Array<T,3> fluidForce, const T& dt, const T& dx){
		Array<T,3> ds = Array<T,3>(0,0,0);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			T toLattice = dt/dx;
			time.push_back(timeLB*dt);
			Array<T,3> f = Array<T,3>(0,0,0);
			Array<T,3> a = Array<T,3>(0,0,0);
			Array<T,3> v = Array<T,3>(0,0,0);
			Array<T,3> v_prev = Array<T,3>(0,0,0);
			Array<T,3> c = Array<T,3>(0,0,0);
			Array<T,3> c_prev = Array<T,3>(0,0,0);
			if(velocity.size()>0){v_prev = velocity.back();}
			if(location.size()>0){c_prev = location.back();}
			for(int i = 0; i<3; i++){
				if(i==2){f[i] += fluidForce[i] + g*mass;}
				else{f[i] += fluidForce[i];}
				a[i] += f[i] / mass;
				v[i] += v_prev[i] + a[i]*dt;
				c[i] += c_prev[i] + v[i]*dt;
				ds[i] += c[i] - c_prev[i];
				ds[i] = ds[i]*toLattice;
			}
			#ifdef PLB_DEBUG
				pcout << "Force on object= "<< array_string(f) <<std::endl;
				pcout << "Acceleration on object= "<< array_string(a) <<std::endl;
				pcout << "Object velocity= "<< array_string(v) <<std::endl;
				pcout << "Object location= "<< array_string(c) <<std::endl;
				mesg = "[DEBUG] DONE Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			force.push_back(f);
			acceleration.push_back(a);
			velocity.push_back(v);
			location.push_back(c);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return ds;
	}



} // NAMESPACE PLB

#endif // VELOCITY_HH
