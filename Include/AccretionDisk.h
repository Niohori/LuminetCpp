#pragma once
#ifndef ACCRETIONDISK_H
#define ACCRETIONDISK_H
#include <algorithm>
#include <cmath>
#include <exception>
#include <iostream>
#include <limits>
#include <memory>
#include <utility>
#include <vector>
#include <map>
#include <list>
#include <iomanip>
#include <unordered_map>
#include <map>
#include <stack>
#include <queue>
#include <deque>
#include <chrono>
#include <random>
#include <ctime>

#include <algorithm>
#include <iostream>
#include <limits>
#include <set>
#include <stdexcept>
#include <vector>
#include <dlib/threads.h>
#include "BlackHolePhysics.h"
#include <memory>


using namespace dlib;
;

namespace accretiondisk{
	// Structure representing an emmitting particle
	//r,phi are the polarcoordinates on the disk and velocity the velocity
		//x,y,  Fs,Fo,redShift are the parameters seen by the observer
	struct Particle {
		double r,phi,x, y, Fs,Fo,redShift, velocity;
		Particle() {
			;
		};
		Particle(double r_, double phi_,double x_P, double y_P, double Fs_P,double Fo_P, double redShift_P, double velocity_) :
			r(r_), 
			phi(phi_),
			x(x_P),
			y(y_P),
			Fs(Fs_P),
			Fo(Fo_P),
			redShift(redShift_P),		
			velocity(velocity_) {
			;
		}
		

	};

	class AccretionDisk {
	public:
		AccretionDisk(void);
		AccretionDisk(const double&, const double&, const double&, const double&, const int&);
		std::vector<Particle> getPrimaryImage();
		std::vector<Particle> getSecundaryImage();
		double getMaxFluxPrimary();
		double getMinFluxPrimary();
		double getMaxFluxSecundary();
		double getMinFluxSecundary();
		void updateDisk(const double&);

	public://variables

	private://methods
		void createDisk();
		double reduceAngle(const double& );


	private://variables
		double M;
		double inclination;
		double innerRadius;
		double outerRadius;
		std::vector<Particle> disk;
		int nParticles;
		int nRadii=100;
		solver_params solver_params_;
		double accretionRate =  1e-8;
		std::vector<Particle> primaryParticles;
		std::vector<Particle> secundaryParticles;
		double maxFluxPrimary; 
		double minFluxPrimary;
		double maxFluxSecundary;
		double minFluxSecundary;
		



	};

} //namespace meshes
#endif