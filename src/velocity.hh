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
	void SurfaceVelocity<T>::initialize(const T& mass_, const T& g_, const T& rho_, const T& dt, const T& dx,
		const ConnectedTriangleSet<T>& triangleSet)
	{
		#ifdef PLB_DEBUG
			std::string mesg = "[DEBUG] Initializing SurfaceVelocity";
			if(master){std::cout << mesg << std::endl;}
			global::log(mesg);
		#endif

		T m3 = dx*dx*dx;
		T t2 = util::sqr(dt);
		mass = mass_;
		rho = rho_ * m3;
		g = -g_*util::sqr(dt)/dx;

		Array<T,3> a = Array<T,3>(0, 0, g);
		acceleration.push_back(a);
		Array<T,3> f = Array<T,3>(0,0,g*mass);
		forceList.push_back(f);
		Array<T,3> m = Array<T,3>(0,0,0);
		torque.push_back(m);

		plint size = triangleSet.getNumVertices();

		velocity.resize(size);
		velocity.reserve(size);

		for(int i = 0; i < size; i++){
			Array<T,3> v = Array<T,3>(0,0,0);
			velocity.push_back(v);
		}

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
	Array<T,6> SurfaceVelocity<T>::gettorqueOfInertia(const Array<T,3>& cg, const ConnectedTriangleSet<T>& triangleSet)
	{
		// Compute the torque of Inertia
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
	Array<T,3> SurfaceVelocity<T>::update(const T& timeLB, const Array<T,3>& force, const Array<T,3>& torque,
		ConnectedTriangleSet<T>& triangleSet)
	{
		Array<T,3> cg = Array<T,3>(0,0,0);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif

			plint n = triangleSet.getNumVertices();
			std::vector<Array<T,3> > oldVertices;
			oldVertices.resize(n);
			oldVertices.reserve(n);

			for(plint i = 0; i<n; i++){
				oldVertices[i] = triangleSet.getVertex(i);
			}

			cg = getCG(oldVertices);

			Array<T,3> f = force;
			f[2] += mass * g;

			Array<T,6> I = gettorqueOfInertia(cg, triangleSet);

			Array<T,3> alpha = getAlpha(torque, I);

			Array<T,3> a = f / mass;

			T dt = timeLB - time.back();
			time.push_back(timeLB);

			Array<T,3> v = a * dt;

			Array<T,3> omega = alpha * dt;

			std::vector<Array<T,3> > newVertices;
			newVertices.resize(n);
			newVertices.reserve(n);

			for(plint i = 0; i < n; i++){
				Array<T,3> u = Array<T,3>(1,1,1);
				Array<T,3> newPosition = getRotatedPosition(oldVertices[i], omega, u, cg);
				verticesVelocity[i] = getRotationalVelocity(oldVertices[i], omega, u, cg) + v;
				for(int r = 0; r < 3; r++){
					newPosition[r] += v[r] * dt;
				}
				newVertices[i] = newPosition;
			}

			triangleSet.swapGeometry(newVertices);

			cg = getCG(newVertices);

			#ifdef PLB_DEBUG
				pcout << "Fluid Force on object = "<< array_string(force) <<std::endl;
				pcout << "torque on object= "<< array_string(torque) <<std::endl;
				pcout << "Acceleration on object= "<< array_string(a) <<std::endl;
				pcout << "Rotational Acceleration on object= "<< array_string(alpha) <<std::endl;
				pcout << "Object velocity= "<< array_string(v) <<std::endl;
				pcout << "Object rotational velocity= "<< array_string(omega) <<std::endl;
				pcout << "Object location= "<< array_string(cg) <<std::endl;
				mesg = "[DEBUG] DONE Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
			#endif
			forceList.push_back(f);
			acceleration.push_back(a);
			velocity.push_back(v);
			location.push_back(cg);
			angular_acceleration.push_back(alpha);
			angular_velocity.push_back(omega);
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return cg;
	}



} // NAMESPACE PLB

#endif // VELOCITY_HH
