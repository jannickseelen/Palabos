#ifndef VELOCITY_HH
#define VELOCITY_HH

#include <cmath>

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

		previous.a_lb = Array<T,3>(0,0,0);
		previous.v_lb = Array<T,3>(0,0,0);
		previous.alpha_lb = Array<T,3>(0,0,0);
		previous.omega_lb = Array<T,3>(0,0,0);

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
				if(iVolume<0){iVolume *= -1; }
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
			#ifdef PLB_DEBUG
				pcout << "[DEBUG]: Moment of Inertia =" <<std::endl;
				pcout << "| "<< Ixx <<" "<<Ixy<<" "<<Ixz<<" |"<<std::endl;
				pcout << "| "<< Ixy <<" "<<Iyy<<" "<<Ixy<<" |"<<std::endl;
				pcout << "| "<< Ixz <<" "<<Ixy<<" "<<Izz<<" |"<<std::endl;
			#endif
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
			if(M[0] != 0 || M[1] != 0 || M[2] != 0){
				T Ixx = I[0];
				T Iyy = I[1];
				T Izz = I[2];
				T Ixy = I[3]; // Ixy == Iyx
				T Ixz = I[4]; // Ixz == Izx
				T Iyz = I[5]; // Iyz == Izy
				T Tx = M[0];
				T Ty = M[1];
				T Tz = M[2];
				alpha[0] = Tx / Ixx + Ty / Ixy + Tz / Ixz;
				alpha[1] = Tx / Ixy + Ty / Iyy + Tz / Iyz;
				alpha[2] = Tx / Ixz + Ty / Iyz + Tz / Izz;
			}
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return alpha;
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::getRotation(const Array<T,3>& vertex, const Array<T,3>& cg, const Array<T,3>& dtheta)
	{
		Array<T,3> newVertex = Array<T,3>(0,0,0);
		try{
			if(dtheta[0] != 0 || dtheta[1] != 0 || dtheta[2] != 0){
				// X, Y, Z Rotation
				const T PI = std::acos(-1);
				// 1. Rotation about the X-axis
				T dz = vertex[2] - cg[2];
				T dy = vertex[1] - cg[1];
				T theta = std::atan2(dz,dy);
				if(theta<0){ theta += 2*PI; }
				T r = std::sqrt(dz*dz + dy*dy);
				T angle = theta + dtheta[0];
				dy += r * std::cos(angle);
				dz += r * std::sin(angle);
				// 2. Rotation about the Y-axis
				T dx = vertex[0] - cg[0];
				theta = std::atan2(dz,dx);
				if(theta<0){ theta += 2*PI; }
				r = std::sqrt(dz*dz + dx*dx);
				angle = theta + dtheta[1];
				dx += r * std::cos(angle);
				dz += r * std::sin(angle);
				// Rotation about the Z-axis
				theta = std::atan2(dx,dy);
				if(theta<0){ theta += 2*PI; }
				angle = theta + dtheta[2];
				r = std::sqrt(dy*dy + dx*dx);
				dy += r * std::cos(angle);
				dx += r * std::sin(angle);
				// Then the new Vertex Becomes
				newVertex[0] = cg[0] + dx;
				newVertex[1] = cg[1] + dy;
				newVertex[2] = cg[2] + dz;
			}
			else{return vertex; }
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return newVertex;
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::getTotalVelocity(const Array<T,3>& vertex, const Array<T,3>& cg, const Array<T,3>& omega_lb,
		const Array<T,3>& v_lb)
	{
		Array<T,3> v = Array<T,3>(0,0,0);
		try{
			if(omega_lb[0] != 0 || omega_lb[1] != 0 || omega_lb[2] !=0){
				// X, Y, Z Rotational Velocity to Linear
				T dx = vertex[0] - cg[0];
				T dy = vertex[1] - cg[1];
				T dz = vertex[2] - cg[2];
				// 1. Angular Velocity about the X-axis
				T r = std::sqrt(dx*dx + dy*dy);
				v[1] += r * std::cos(omega_lb[0]);
				v[2] += r * std::sin(omega_lb[0]);
				// 2. Rotation about the Y-axis
				r = std::sqrt(dz*dz + dx*dx);
				v[0] += r * std::cos(omega_lb[1]);
				v[2] += r * std::sin(omega_lb[1]);
				// Rotation about the Z-axis
				r = std::sqrt(dy*dy + dx*dx);
				v[1] += r * std::cos(omega_lb[2]);
				v[0] += r * std::sin(omega_lb[2]);
				// Then the new Vertex Becomes
				v += v_lb;
			}
			else{ return v_lb; }
		}
		catch(const std::exception& e){exHandler(e,__FILE__,__FUNCTION__,__LINE__);}
		return v;
	}

	template<typename T>
	Array<T,3> SurfaceVelocity<T>::update(const IncomprFlowParam<T>& p, const T& time_lb, const Array<T,3>& force,
		const Array<T,3>& torque, ConnectedTriangleSet<T>& triangleSet)
	{
		Array<T,3> cg_lb = Array<T,3>(0,0,0);
		try{
			#ifdef PLB_DEBUG
				std::string mesg = "[DEBUG] Updating SurfaceVelocity";
				if(master){std::cout << mesg << std::endl;}
				global::log(mesg);
				pcout << "Input in Dimensionless Units" << std::endl;
				pcout << "FluidForce= "<< array_string(force) << std::endl;
				pcout << "FluidTorque= "<< array_string(torque) << std::endl;
			#endif

			const T dt = p.getDeltaT();
			const T dx = p.getDeltaX();

			plint n = triangleSet.getNumVertices();
			std::vector<Array<T,3> > oldVertices;
			oldVertices.resize(n);
			oldVertices.reserve(n);

			for(plint i = 0; i<n; i++){oldVertices[i] = triangleSet.getVertex(i);}

			cg_lb = getCG(oldVertices);

			T mass_lb = mass / (dx*dx*dx);
			T g_lb = g * dt * dt  / dx;
			T gravityForce = mass_lb * g_lb;

			Array<T,3> f_lb = Array<T,3>(0,0,0);
			f_lb = force; // * dt*dt / (dx * dx * dx *dx);
			f_lb[2] += gravityForce;

			Array<T,3> torque_lb = Array<T,3>(0,0,0);
			torque_lb = torque; // * dt*dt / (dx * dx * dx * dx * dx);

			Array<T,6> I_lb = Array<T,6>(0,0,0,0,0,0);
			I_lb = getMomentOfInertia(cg_lb, triangleSet);

			Array<T,3> a_lb = Array<T,3>(0,0,0);
			a_lb = previous.a_lb + f_lb / mass_lb;

			Array<T,3> v_lb = Array<T,3>(0,0,0);
			v_lb = previous.v_lb + a_lb * (T)1.0;

			Array<T,3> ds_lb = Array<T,3>(0,0,0);
			ds_lb = previous.v_lb * (T)1.0 + (T)0.5 * a_lb * (T)1.0 * (T)1.0;

			pcout << "ds LB= " << array_string(ds_lb) << std::endl;

			Array<T,3> alpha_lb = Array<T,3>(0,0,0);
			alpha_lb = previous.alpha_lb + getAlpha(torque_lb, I_lb);

			Array<T,3> omega_lb = Array<T,3>(0,0,0);
			omega_lb = previous.omega_lb + alpha_lb * (T)1.0;

			Array<T,3> dtheta_lb = Array<T,3>(0,0,0);
			dtheta_lb = previous.omega_lb * (T)1.0 + (T)0.5 * alpha_lb * (T)1.0 * (T)1.0;

			std::vector<Array<T,3> > newVertices;
			newVertices.resize(n);
			newVertices.reserve(n);

			verticesVelocity.clear();
			verticesVelocity.resize(n);
			verticesVelocity.reserve(n);

			for(plint i = 0; i < n; i++){
				//pcout << "old Vertex= " << array_string(oldVertices[i]);
				newVertices[i] = getRotation(oldVertices[i],cg_lb,dtheta_lb);
				newVertices[i] += ds_lb;
				//pcout << " new Vertex= " << array_string(newVertices[i]);
				verticesVelocity[i] = getTotalVelocity(oldVertices[i],cg_lb,omega_lb,v_lb);
				//pcout << " Vertex velocity=  "<< array_string(verticesVelocity[i]) << std::endl;
			}

			triangleSet.swapGeometry(newVertices);

			cg_lb = getCG(newVertices);

			previous.v_lb = v_lb;
			previous.a_lb = a_lb;
			previous.omega_lb = omega_lb;
			previous.alpha_lb = alpha_lb;

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
				pcout << "Object location in lattice units= "<< array_string(cg_lb) <<std::endl;
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
