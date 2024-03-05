#pragma once
#ifndef ISOREDSHIFTS_H
#define ISOREDSHIFTS_H
#include <iostream>
#include <cmath>
#include <fstream>
#include <random>
#include <algorithm>
#include <unordered_map>
#include <iostream>
#include <cmath>
#include <vector>
#include <cmath>
#include <functional>
#include <map>
#include <numeric>
#include <iomanip>

#include "TensorCalculus.h"
#include "BlackHolePhysics.h"
#include "IsoRadials.h"
#include "utilities.h"
#include "mesh.h"
#include "isolines.h"

#include <dlib/threads.h>

using namespace dlib;
using namespace meshes;

#define xsect(p1,p2) (h[p2]*xh[p1]-h[p1]*xh[p2])/(h[p2]-h[p1])
#define ysect(p1,p2) (h[p2]*yh[p1]-h[p1]*yh[p2])/(h[p2]-h[p1])

class IsoRedShift {
public://methods
	IsoRedShift();
	IsoRedShift(const double& angle, const double& bh_mass, const double& lower_radius, const double& upper_radius, const size_t& _n_radii_, const size_t& _n_angles_, const double&);
	IsoRedShift(const double&, const double&, const double&, const std::map<double, std::pair<int, Isoradial*> >&);
	~IsoRedShift();
	std::multimap<double, std::vector<meshes::Point> > get_isolines(const size_t n);
	std::pair<std::vector<double>, std::vector<double>> get_ISCO_curve();
	std::pair<std::vector<double>, std::vector<double>> get_ConcaveHull();
	double  getMaxX();
	
private://variables
	double theta_0;  // Inclination
	double redshift;
	double M;  // Black hole mass
	double rS;
	double rIsco;
	BHphysics BHp;
	std::vector<double> angles;
	std::vector<double> radii;
	size_t n_radii;
	size_t n_angles;
	double lower_radius;
	double upper_radius;
	std::vector<double> x;
	std::vector<double> y;
	double max_radius;
	irs_solver_params irs_solver_params_;
	find_redshift_params find_redshift_params_;
	angular_properties angular_properties_;
	solver_params solver_params_;
	plot_params plot_params_;
	ir_params ir_params_;
	std::pair<std::vector<double>, std::vector<double> > ISCO_boundary;
	double maxCoordinate;

private://methods
	void make_grid();
	double calculateDistance(const meshes::Point& p1, const meshes::Point& p2);
	double findSmallestDistance(const std::vector<meshes::Point>&);
	

public://variables ? make get method?
	double redshift_;
	double x_;
	double y_;
	double x_max = -100000000000000.0;
	double x_min = -x_max;
	double y_max = x_max;
	double y_min = x_min;
	double redshift_max = -1000000000000.0;
	double redshift_min = -redshift_max;
	std::vector<double> xIsco;
	std::vector<double> yIsco;
	std::vector<meshes::Point> ConcaveHull;
	std::vector<meshes::Point> ISCO;
	double delta = 5.0e-2;
	double redshift_treshold;//TEMPORARY
	std::vector<double> xCoordinates;
	std::vector<double > yCoordinates;
	std::vector<double > redshifts;
	std::vector<double> xGrid;//TEMPORARY
	std::vector<double> yGrid;//TEMPORARY
	std::vector<std::vector<double>> redshiftGrid;//TEMPORARY
};
#endif