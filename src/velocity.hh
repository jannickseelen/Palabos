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
	Array<T,3> SurfaceVelocity<T>::operator()(pluint id)
	{
		return verticesVelocity[id];
	}

	template<typename T>
	void SurfaceVelocity<T>::initialize(const T& mass_, const T& g_, const T& rho_)
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Initializing SurfaceVelocity";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif

		mass = mass_;
		rho = rho_;
		g = -g;

		Array<T,3> a = Array<T,3>(0, 0, g);
		acceleration.push_back(a);
		Array<T,3> f = Array<T,3>(0,0,g*mass);
		forceList.push_back(f);
		Array<T,3> m = Array<T,3>(0,0,0);
		torque.push_back(m);
		T t = 0;
		time.push_back(t);

		#ifdef PLB_DEBUG
			mesg = "[DEBUG] DONE Initializing SurfaceVelocity";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::getCG(std::vector<Array<T,3> > vertexList)
	{
		// Compute the center of gravity
		Array<T,3> cg = Array<T,3>(0,0,0);
		try{
			T x = 0;
			T y = 0;
			T z = 0;
			int size = vertexList.size();

			for(int i = 0; i<size; i++){
				Array<T,3> iVertex = vertexList[i];
				x += iVertex[0];
				y += iVertex[1];
				z += iVertex[2];
			}

			cg = Array<T,3>(x/size, y/size, z/size);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return cg;
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::getArm(const Array<T,3>& p1, const Array<T,3>& p2)
	{
		// Compute the arm r between the position and cg
		Array<T,3> arm = Array<T,3>(0,0,0);
		try{
			for(int i = 0; i<3; i++){
				if(p1[i] > p2[i]){
					if(p1[i] >= 0){
						arm[i] = p1[i] - p2[i];
					}
					if(p1[i] < 0 && p2[i] < 0){
						T a = -1*p2[i];
						T b = -1*p1[i];
						arm[i] = a - b;
					}
				}
				else{
					if(p2[i] >= 0){
						arm[i] = p2[i] - p1[i];
					}
					if(p1[i] < 0 && p2[i] < 0){
						T a = -1*p1[i];
						T b = -1*p2[i];
						arm[i] = a - b;
					}
				}
				if(arm[i] < 0){ arm[i]*= -1; }
			}
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return arm;
	}

	template<typename T>
	Array<T,6> SurfaceVelocity<T>::getMomentOfInertia(const Array<T,3>& cg, const ConnectedTriangleSet<T>& triangleSet)
	{
		// Compute the moment of Inertia
		Array<T,6> I = Array<T,6>(0,0,0,0,0,0);
		try{
			T Ixx = 0;
			T Iyy = 0;
			T Izz = 0;
			T Ixy = 0;
			T Ixz = 0;
			T Iyz = 0;
			T numTriangles = triangleSet.getNumTriangles();
			for(int i = 0; i<numTriangles; i++){
				Array<plint,3> iTriangle = triangleSet.getTriangle(i);
				Array<T,3> a = triangleSet.getVertex(iTriangle[0]);
				Array<T,3> b = triangleSet.getVertex(iTriangle[1]);
				Array<T,3> c = triangleSet.getVertex(iTriangle[2]);
				Array<T,3> d = cg;
				T iVolume = computeTetrahedronSignedVolume(a,b,c,d);
				T iMass = iVolume * rho;
				Array<T,3> iCG = Array<T,3>(0,0,0);
				for(int i = 0; i<3; i++){
					T temp = a[i]+b[i]+c[i]+d[i];
					iCG[i] = temp/4;
				}
				T x = iCG[0] - cg[0];
				T y = iCG[1] - cg[1];
				T z = iCG[2] - cg[2];
				T temp = y*y + z*z;
				Ixx += iMass * temp;
				temp = x*x + z*z;
				Iyy += iMass * temp;
				temp = x*x + y*y;
				Izz += iMass * temp;
				Ixy += iMass * x * y;
				Ixz += iMass * x * z;
				Iyz += iMass * y * z;
			}
			Ixy *= -1;
			Ixz *= -1;
			Iyz *= -1;
			I[0] = Ixx;
			I[1] = Iyy;
			I[2] = Izz;
			I[3] = Ixy;
			I[4] = Ixz;
			I[5] = Iyz;
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return I;
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::getAlpha(const Array<T,3>& M, const Array<T,6>& I)
	{
		// Compute the angular acceleration
		Array<T,3> alpha = Array<T,3>(0,0,0);
		try{
			T x = 0;
			T y = 0;
			T z = 0;
			T denum = I[0]*I[0]*I[1] - I[3]*I[3];
			T temp = I[0]*denum;
			T lhs_a = I[4]*I[4]*I[3]*I[3]/temp;
			T lhs_b = I[3]*I[4]*I[5]/denum;
			T lhs_c = I[4]*I[4]/I[0];
			T lhs_d = -I[3]*I[4]*I[5]/denum;
			T lhs_e = I[0]*I[5]*I[5]/denum;
			T lhs = lhs_a + lhs_b + lhs_c + lhs_d + lhs_e - I[2];
			T rhs_a = I[4]*M[0]/I[0];
			T rhs_b = -I[3]*I[4]*M[1]/denum;
			T rhs_c = I[3]*I[3]*I[4]*M[0]/temp;
			T rhs_d = I[0]*I[5]*M[1]/denum;
			T rhs_e = -I[3]*I[4]*M[0]/denum;
			T rhs = rhs_a + rhs_b + rhs_c + rhs_d + rhs_e - M[2];
			z = rhs/lhs;
			y = I[0]*M[1] - I[3]*M[0] + I[4]*I[5]*z - I[0]*I[5]*z / denum;
			x = M[0] - I[3]*y - I[4]*z / I[0];
			alpha[0] = x;
			alpha[1] = y;
			alpha[2] = z;
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return alpha;
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::update(const IncomprFlowParam<T>& p, const T& time_lb, const Array<T,3>& force_lb,
		const Array<T,3>& torque_lb, ConnectedTriangleSet<T>& triangleSet)
	{
		Array<T,3> cg_lb = Array<T,3>(0,0,0);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				pcout << "Input in Dimensionless Units" << std::endl;
				pcout << "FluidForce= "<< array_string(force_lb) << std::endl;
				pcout << "FluidTorque= "<< array_string(torque_lb) << std::endl;
			#endif

			const T dt = p.getDeltaT();
			const T dx = p.getDeltaX();

			plint n = triangleSet.getNumVertices();
			std::vector<Array<T,3> > oldVertices;
			oldVertices.resize(n);
			oldVertices.reserve(n);

			for(plint i = 0; i<n; i++){
				oldVertices[i] = triangleSet.getVertex(i);
			}

			cg_lb = getCG(oldVertices);

			T mass_lb = mass / (dx*dx*dx);
			T g_lb = g * dt * dt  / dx;
			T gravityForce = mass_lb * g_lb;

			Array<T,3> f_lb = force_lb;
			f_lb[2] += gravityForce;

			Array<T,6> I_lb = getMomentOfInertia(cg_lb, triangleSet);

			Array<T,3> alpha_lb = getAlpha(torque_lb, I_lb);

			Array<T,3> a_lb = f_lb / mass_lb;

			Array<T,3> v_lb = a_lb * (T)1.0;

			Array<T,3> omega_lb = alpha_lb * (T)1.0;

			std::vector<Array<T,3> > newVertices;
			newVertices.resize(n);
			newVertices.reserve(n);
			verticesVelocity.clear();
			verticesVelocity.resize(n);
			verticesVelocity.reserve(n);

			for(plint i = 0; i < n; i++){
				Array<T,3> u = Array<T,3>(1,1,1);
				Array<T,3> newPosition = getRotatedPosition(triangleSet.getVertex(i), omega_lb, u, cg_lb);
				verticesVelocity[i] = getRotationalVelocity(triangleSet.getVertex(i), omega_lb, u, cg_lb) + v_lb;
				for(int r = 0; r < 3; r++){
					newPosition[r] += v_lb[r] * (T)1.0;
				}
				newVertices[i] = newPosition;
			}

			triangleSet.swapGeometry(newVertices);

			cg_lb = getCG(newVertices);

			Array<T,3> f =  f_lb*dx*dx*dx*dx/(dt*dt);
			Array<T,3> t = torque_lb*dx*dx*dx*dx*dx/(dt*dt);
			Array<T,3> v = v_lb*dx/dt;
			Array<T,3> a = a_lb*dx/(dt*dt);
			Array<T,3> alpha = alpha_lb*dx/(dt*dt);
			Array<T,3> omega = omega_lb*dx/dt;
			Array<T,3> cg = cg_lb*dx;

			#ifdef PLB_DEBUG
				pcout << "Kinematics in Physical Units" << std::endl;
				pcout << "Force on object = "<< array_string(f) <<std::endl;
				pcout << "Torque on object= "<< array_string(t) <<std::endl;
				pcout << "Acceleration on object= "<< array_string(a) <<std::endl;
				pcout << "Rotational Acceleration on object= "<< array_string(alpha) <<std::endl;
				pcout << "Object velocity= "<< array_string(v) <<std::endl;
				pcout << "Object rotational velocity= "<< array_string(omega) <<std::endl;
				pcout << "Object location in lattice units= "<< array_string(cg) <<std::endl;
				mesg = "[DEBUG] DONE Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			time.push_back(time_lb);
			forceList.push_back(f);
			acceleration.push_back(a);
			velocity.push_back(v);
			location.push_back(cg);
			angular_acceleration.push_back(alpha);
			angular_velocity.push_back(omega);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return cg_lb;
	}



} // NAMESPACE PLB

#endif // VELOCITY_HH
